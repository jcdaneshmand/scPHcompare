#![forbid(unsafe_op_in_unsafe_fn)]

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::slice;

pub const ENGINE_VERSION_V1: u32 = 1;
pub const STATUS_OK: i32 = 0;
pub const STATUS_NULL_OUTPUT: i32 = 1;
pub const STATUS_NULL_INPUT: i32 = 2;
pub const STATUS_INVALID_INPUT: i32 = 3;
pub const STATUS_PANIC: i32 = 4;
pub const STATUS_NONFINITE_OUTPUT: i32 = 5;
pub const STATUS_UNSUPPORTED_DIMENSION: i32 = 6;

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ScphLandscapeResultV1 {
    pub squared_distance: f64,
    pub active_levels: u64,
    pub event_segments: u64,
    pub first_finite_intervals: u64,
    pub second_finite_intervals: u64,
    pub engine_version: u32,
    pub status: i32,
}

impl ScphLandscapeResultV1 {
    fn with_status(status: i32) -> Self {
        Self {
            squared_distance: f64::NAN,
            active_levels: 0,
            event_segments: 0,
            first_finite_intervals: 0,
            second_finite_intervals: 0,
            engine_version: ENGINE_VERSION_V1,
            status,
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Point {
    x: f64,
    y: f64,
}

#[derive(Clone, Copy, Debug)]
struct OrderedF64(f64);

impl PartialEq for OrderedF64 {
    fn eq(&self, other: &Self) -> bool {
        self.0.total_cmp(&other.0) == Ordering::Equal
    }
}

impl Eq for OrderedF64 {}

impl PartialOrd for OrderedF64 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OrderedF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

#[derive(Debug)]
struct MaxTree {
    leaf_count: usize,
    values: Vec<f64>,
}

impl MaxTree {
    fn new(count: usize) -> Self {
        let leaf_count = count.max(1).next_power_of_two();
        Self {
            leaf_count,
            values: vec![f64::NEG_INFINITY; 2 * leaf_count],
        }
    }

    fn update(&mut self, index: usize, value: f64) {
        let mut node = self.leaf_count + index;
        self.values[node] = value;
        while node > 1 {
            node /= 2;
            self.values[node] = self.values[2 * node].max(self.values[2 * node + 1]);
        }
    }

    fn first_above(&self, threshold: f64) -> Option<usize> {
        if self.values[1] <= threshold {
            return None;
        }
        let mut node = 1;
        while node < self.leaf_count {
            let left = 2 * node;
            node = if self.values[left] > threshold {
                left
            } else {
                left + 1
            };
        }
        Some(node - self.leaf_count)
    }
}

#[derive(Debug)]
struct IntervalPool {
    births: Vec<f64>,
    deaths: Vec<BTreeMap<OrderedF64, usize>>,
    maxima: MaxTree,
    remaining: usize,
}

impl IntervalPool {
    fn new(intervals: &[(f64, f64)]) -> Self {
        let mut births: Vec<f64> = intervals.iter().map(|item| item.0).collect();
        births.sort_by(f64::total_cmp);
        births.dedup_by(|a, b| a.total_cmp(b) == Ordering::Equal);
        let mut deaths = vec![BTreeMap::new(); births.len()];
        for &(birth, death) in intervals {
            let index = births
                .binary_search_by(|value| value.total_cmp(&birth))
                .expect("birth originated from the indexed interval set");
            *deaths[index].entry(OrderedF64(death)).or_insert(0) += 1;
        }
        let mut maxima = MaxTree::new(births.len());
        for (index, group) in deaths.iter().enumerate() {
            let maximum = group
                .last_key_value()
                .map_or(f64::NEG_INFINITY, |(value, _)| value.0);
            maxima.update(index, maximum);
        }
        Self {
            births,
            deaths,
            maxima,
            remaining: intervals.len(),
        }
    }

    fn refresh(&mut self, index: usize) {
        let maximum = self.deaths[index]
            .last_key_value()
            .map_or(f64::NEG_INFINITY, |(value, _)| value.0);
        self.maxima.update(index, maximum);
    }

    fn pop_first(&mut self) -> Option<(f64, f64)> {
        if self.remaining == 0 {
            return None;
        }
        let index = self.maxima.first_above(f64::NEG_INFINITY)?;
        let key = *self.deaths[index].last_key_value()?.0;
        let count = self.deaths[index].get_mut(&key)?;
        *count -= 1;
        if *count == 0 {
            self.deaths[index].remove(&key);
        }
        self.remaining -= 1;
        self.refresh(index);
        Some((self.births[index], key.0))
    }

    fn pop_first_with_death_above(&mut self, threshold: f64) -> Option<(f64, f64)> {
        let index = self.maxima.first_above(threshold)?;
        let key = *self.deaths[index].last_key_value()?.0;
        let count = self.deaths[index].get_mut(&key)?;
        *count -= 1;
        if *count == 0 {
            self.deaths[index].remove(&key);
        }
        self.remaining -= 1;
        self.refresh(index);
        Some((self.births[index], key.0))
    }

    fn insert(&mut self, birth: f64, death: f64) {
        let index = self
            .births
            .binary_search_by(|value| value.total_cmp(&birth))
            .expect("residual birth originated from the indexed interval set");
        *self.deaths[index].entry(OrderedF64(death)).or_insert(0) += 1;
        self.remaining += 1;
        self.refresh(index);
    }
}

fn critical_landscape(intervals: &[(f64, f64)]) -> Vec<Vec<Point>> {
    let mut pool = IntervalPool::new(intervals);
    let mut landscape = Vec::new();
    // Residual intervals created below may be identical even when their source
    // barcode intervals were not.  Collapsing those residual duplicates and
    // cloning the completed suffix is invalid: the residuals can belong to
    // different earlier landscape chains.  Consume exactly one interval per
    // outer level so every chain is constructed independently.
    while let Some((birth, mut death)) = pool.pop_first() {
        let mut level = Vec::new();
        level.push(Point { x: birth, y: 0.0 });
        level.push(Point {
            x: (birth + death) / 2.0,
            y: (death - birth) / 2.0,
        });
        loop {
            let Some((next_birth, next_death)) = pool.pop_first_with_death_above(death) else {
                level.push(Point { x: death, y: 0.0 });
                break;
            };
            if next_birth > death {
                level.push(Point { x: death, y: 0.0 });
            }
            if next_birth >= death {
                level.push(Point {
                    x: next_birth,
                    y: 0.0,
                });
            } else {
                level.push(Point {
                    x: (next_birth + death) / 2.0,
                    y: (death - next_birth) / 2.0,
                });
                pool.insert(next_birth, death);
            }
            level.push(Point {
                x: (next_birth + next_death) / 2.0,
                y: (next_death - next_birth) / 2.0,
            });
            death = next_death;
        }
        landscape.push(level);
    }
    landscape
}

#[derive(Debug, Default)]
struct KahanSum {
    sum: f64,
    correction: f64,
}

impl KahanSum {
    fn add(&mut self, value: f64) {
        let adjusted = value - self.correction;
        let updated = self.sum + adjusted;
        self.correction = (updated - self.sum) - adjusted;
        self.sum = updated;
    }
}

fn value_at(points: &[Point], location: f64, segment: &mut usize) -> f64 {
    if points.is_empty() || location < points[0].x || location > points[points.len() - 1].x {
        return 0.0;
    }
    while *segment + 1 < points.len() && points[*segment + 1].x <= location {
        *segment += 1;
    }
    if *segment + 1 == points.len() || points[*segment].x == location {
        return points[*segment].y;
    }
    let left = points[*segment];
    let right = points[*segment + 1];
    if right.x == left.x {
        return right.y;
    }
    left.y + (right.y - left.y) * (location - left.x) / (right.x - left.x)
}

fn integrate_level(first: &[Point], second: &[Point], total: &mut KahanSum) -> u64 {
    let mut first_node = 0;
    let mut second_node = 0;
    let mut first_segment = 0;
    let mut second_segment = 0;
    let mut previous_x: Option<f64> = None;
    let mut previous_difference = 0.0;
    let mut segments = 0;
    while first_node < first.len() || second_node < second.len() {
        let location = match (first.get(first_node), second.get(second_node)) {
            (Some(a), Some(b)) => match a.x.total_cmp(&b.x) {
                Ordering::Less => {
                    first_node += 1;
                    a.x
                }
                Ordering::Greater => {
                    second_node += 1;
                    b.x
                }
                Ordering::Equal => {
                    first_node += 1;
                    second_node += 1;
                    a.x
                }
            },
            (Some(a), None) => {
                first_node += 1;
                a.x
            }
            (None, Some(b)) => {
                second_node += 1;
                b.x
            }
            (None, None) => break,
        };
        let difference = value_at(first, location, &mut first_segment)
            - value_at(second, location, &mut second_segment);
        if let Some(left) = previous_x {
            let width = location - left;
            if width > 0.0 {
                total.add(
                    width
                        * (previous_difference * previous_difference
                            + previous_difference * difference
                            + difference * difference)
                        / 3.0,
                );
                segments += 1;
            }
        }
        previous_x = Some(location);
        previous_difference = difference;
    }
    segments
}

fn validate_intervals(births: &[f64], deaths: &[f64]) -> Result<Vec<(f64, f64)>, i32> {
    if births.len() != deaths.len() {
        return Err(STATUS_INVALID_INPUT);
    }
    let mut intervals = Vec::with_capacity(births.len());
    for (&birth, &death) in births.iter().zip(deaths) {
        if !birth.is_finite() || !death.is_finite() || death <= birth {
            return Err(STATUS_INVALID_INPUT);
        }
        intervals.push((birth, death));
    }
    Ok(intervals)
}

pub fn landscape_squared_l2(
    first_births: &[f64],
    first_deaths: &[f64],
    second_births: &[f64],
    second_deaths: &[f64],
    dimension: u32,
) -> Result<ScphLandscapeResultV1, i32> {
    if dimension > 1 {
        return Err(STATUS_UNSUPPORTED_DIMENSION);
    }
    let first = validate_intervals(first_births, first_deaths)?;
    let second = validate_intervals(second_births, second_deaths)?;
    let first_landscape = critical_landscape(&first);
    let second_landscape = critical_landscape(&second);
    let active_levels = first_landscape.len().max(second_landscape.len());
    let mut squared = KahanSum::default();
    let mut event_segments = 0;
    for level in 0..active_levels {
        let a = first_landscape.get(level).map_or(&[][..], Vec::as_slice);
        let b = second_landscape.get(level).map_or(&[][..], Vec::as_slice);
        event_segments += integrate_level(a, b, &mut squared);
    }
    if !squared.sum.is_finite() || squared.sum < -1e-12 {
        return Err(STATUS_NONFINITE_OUTPUT);
    }
    Ok(ScphLandscapeResultV1 {
        squared_distance: squared.sum.max(0.0),
        active_levels: active_levels as u64,
        event_segments,
        first_finite_intervals: first.len() as u64,
        second_finite_intervals: second.len() as u64,
        engine_version: ENGINE_VERSION_V1,
        status: STATUS_OK,
    })
}

#[unsafe(no_mangle)]
#[allow(unsafe_code)]
/// Compute one exact persistence-landscape squared-L2 distance through C ABI v1.
///
/// # Safety
///
/// `output` must point to writable `ScphLandscapeResultV1` storage. Each input
/// pointer with a nonzero corresponding length must reference that many valid,
/// contiguous `f64` values for the duration of the call. Input buffers are read
/// only and must not overlap writable output storage.
pub unsafe extern "C" fn scph_landscape_l2_v1(
    first_births: *const f64,
    first_deaths: *const f64,
    first_len: usize,
    second_births: *const f64,
    second_deaths: *const f64,
    second_len: usize,
    dimension: u32,
    output: *mut ScphLandscapeResultV1,
) -> i32 {
    if output.is_null() {
        return STATUS_NULL_OUTPUT;
    }
    if (first_len > 0 && (first_births.is_null() || first_deaths.is_null()))
        || (second_len > 0 && (second_births.is_null() || second_deaths.is_null()))
    {
        unsafe { output.write(ScphLandscapeResultV1::with_status(STATUS_NULL_INPUT)) };
        return STATUS_NULL_INPUT;
    }
    let first_birth_slice = if first_len == 0 {
        &[]
    } else {
        unsafe { slice::from_raw_parts(first_births, first_len) }
    };
    let first_death_slice = if first_len == 0 {
        &[]
    } else {
        unsafe { slice::from_raw_parts(first_deaths, first_len) }
    };
    let second_birth_slice = if second_len == 0 {
        &[]
    } else {
        unsafe { slice::from_raw_parts(second_births, second_len) }
    };
    let second_death_slice = if second_len == 0 {
        &[]
    } else {
        unsafe { slice::from_raw_parts(second_deaths, second_len) }
    };
    let outcome = catch_unwind(AssertUnwindSafe(|| {
        landscape_squared_l2(
            first_birth_slice,
            first_death_slice,
            second_birth_slice,
            second_death_slice,
            dimension,
        )
    }));
    let result = match outcome {
        Ok(Ok(result)) => result,
        Ok(Err(status)) => ScphLandscapeResultV1::with_status(status),
        Err(_) => ScphLandscapeResultV1::with_status(STATUS_PANIC),
    };
    let status = result.status;
    unsafe { output.write(result) };
    status
}

#[unsafe(no_mangle)]
#[allow(unsafe_code, clippy::too_many_arguments)]
/// R `.C` compatibility shim for the versioned C ABI.
///
/// # Safety
///
/// All scalar pointers must reference writable or readable storage matching
/// their declared types. Each interval pointer with a positive scalar length
/// must reference that many contiguous `f64` values. R owns every buffer and
/// must keep it alive for the duration of the call.
pub unsafe extern "C" fn scph_landscape_l2_r_v1(
    first_births: *const f64,
    first_deaths: *const f64,
    first_len: *const i32,
    second_births: *const f64,
    second_deaths: *const f64,
    second_len: *const i32,
    dimension: *const i32,
    squared_distance: *mut f64,
    active_levels: *mut f64,
    event_segments: *mut f64,
    first_finite_intervals: *mut f64,
    second_finite_intervals: *mut f64,
    engine_version: *mut i32,
    status: *mut i32,
) {
    if first_len.is_null()
        || second_len.is_null()
        || dimension.is_null()
        || squared_distance.is_null()
        || active_levels.is_null()
        || event_segments.is_null()
        || first_finite_intervals.is_null()
        || second_finite_intervals.is_null()
        || engine_version.is_null()
        || status.is_null()
    {
        return;
    }
    let first_count = unsafe { *first_len };
    let second_count = unsafe { *second_len };
    let requested_dimension = unsafe { *dimension };
    let mut result = ScphLandscapeResultV1::with_status(STATUS_INVALID_INPUT);
    if first_count >= 0 && second_count >= 0 && requested_dimension >= 0 {
        unsafe {
            scph_landscape_l2_v1(
                first_births,
                first_deaths,
                first_count as usize,
                second_births,
                second_deaths,
                second_count as usize,
                requested_dimension as u32,
                &mut result,
            );
        }
    }
    unsafe {
        squared_distance.write(result.squared_distance);
        active_levels.write(result.active_levels as f64);
        event_segments.write(result.event_segments as f64);
        first_finite_intervals.write(result.first_finite_intervals as f64);
        second_finite_intervals.write(result.second_finite_intervals as f64);
        engine_version.write(result.engine_version as i32);
        status.write(result.status);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn squared(first: &[(f64, f64)], second: &[(f64, f64)]) -> f64 {
        let first_births: Vec<_> = first.iter().map(|value| value.0).collect();
        let first_deaths: Vec<_> = first.iter().map(|value| value.1).collect();
        let second_births: Vec<_> = second.iter().map(|value| value.0).collect();
        let second_deaths: Vec<_> = second.iter().map(|value| value.1).collect();
        landscape_squared_l2(
            &first_births,
            &first_deaths,
            &second_births,
            &second_deaths,
            0,
        )
        .unwrap()
        .squared_distance
    }

    #[test]
    fn analytical_fixtures_are_exact() {
        assert!((squared(&[(0.0, 2.0)], &[]) - 2.0 / 3.0).abs() <= 1e-15);
        assert!((squared(&[(0.0, 2.0)], &[(0.25, 2.25)]) - 7.0 / 64.0).abs() <= 1e-15);
        assert!((squared(&[(0.499, 0.501)], &[]) - 0.002_f64.powi(3) / 12.0).abs() <= 1e-20);
    }

    #[test]
    fn duplicates_create_consecutive_levels() {
        assert!((squared(&[(0.0, 2.0), (0.0, 2.0)], &[]) - 4.0 / 3.0).abs() <= 1e-15);
    }

    #[test]
    fn residual_duplicates_do_not_clone_incomplete_levels() {
        let intervals = [(1.0, 4.0), (3.0, 4.0), (0.0, 2.0), (1.0, 2.0)];
        let expected_norm: f64 = intervals
            .iter()
            .map(|(birth, death)| (death - birth) * (death - birth) * (death - birth) / 12.0)
            .sum();
        let result = landscape_squared_l2(
            &intervals.iter().map(|item| item.0).collect::<Vec<_>>(),
            &intervals.iter().map(|item| item.1).collect::<Vec<_>>(),
            &[],
            &[],
            1,
        )
        .unwrap();
        assert_eq!(result.active_levels, 3);
        assert!((result.squared_distance - expected_norm).abs() <= 1e-15);
    }

    #[test]
    fn all_small_four_interval_norms_and_depths_are_exact() {
        let mut universe = Vec::new();
        for birth in 0..6 {
            for death in (birth + 1)..6 {
                universe.push((birth as f64, death as f64));
            }
        }
        for first in 0..universe.len() {
            for second in (first + 1)..universe.len() {
                for third in (second + 1)..universe.len() {
                    for fourth in (third + 1)..universe.len() {
                        let intervals = [
                            universe[first],
                            universe[second],
                            universe[third],
                            universe[fourth],
                        ];
                        let births: Vec<_> = intervals.iter().map(|item| item.0).collect();
                        let deaths: Vec<_> = intervals.iter().map(|item| item.1).collect();
                        let result = landscape_squared_l2(&births, &deaths, &[], &[], 1).unwrap();
                        let expected_norm: f64 = intervals
                            .iter()
                            .map(|(birth, death)| {
                                (death - birth) * (death - birth) * (death - birth) / 12.0
                            })
                            .sum();
                        let expected_depth = (0..5)
                            .map(|segment| {
                                let location = segment as f64 + 0.5;
                                intervals
                                    .iter()
                                    .filter(|(birth, death)| *birth < location && location < *death)
                                    .count()
                            })
                            .max()
                            .unwrap();
                        assert_eq!(result.active_levels, expected_depth as u64);
                        assert!((result.squared_distance - expected_norm).abs() <= 1e-12);
                    }
                }
            }
        }
    }

    #[test]
    fn result_is_symmetric_and_dimension_separated() {
        let a = [(0.0, 2.0), (0.5, 3.0), (1.0, 1.5)];
        let b = [(0.25, 2.25), (0.75, 2.5)];
        assert_eq!(squared(&a, &b).to_bits(), squared(&b, &a).to_bits());
        assert_eq!(
            landscape_squared_l2(&[0.0], &[2.0], &[], &[], 2).unwrap_err(),
            STATUS_UNSUPPORTED_DIMENSION
        );
    }

    #[test]
    fn invalid_intervals_are_rejected() {
        assert_eq!(
            landscape_squared_l2(&[0.0], &[0.0], &[], &[], 0).unwrap_err(),
            STATUS_INVALID_INPUT
        );
        assert_eq!(
            landscape_squared_l2(&[f64::NAN], &[1.0], &[], &[], 0).unwrap_err(),
            STATUS_INVALID_INPUT
        );
    }
}
