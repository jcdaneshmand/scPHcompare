test_that("MV17-G full calibration seeds and grouped queues are complete", {
  seeds<-mv17g_seed_registry_v1();expect_equal(nrow(seeds),91476L);expect_equal(length(unique(seeds$seed)),91476L);expect_equal(range(seeds$seed),c(174001L,265476L));expect_equal(as.integer(table(seeds$view)),c(39204L,52272L));primary<-mv17g_grouped_queue_v1("primary");expect_equal(nrow(primary),1188L);expect_equal(sum(primary$scientific_runs),91740L);expect_equal(sum(primary$null_family=="observed"),264L);expect_equal(sum(primary$replicate_count==99L),924L)
})

test_that("MV17-G repeat is a seed-exact six-unit subset", {
  repeat_orders<-data.frame(view=rep(c("cell","gene"),each=3L),unit_order=rep(c(1L,66L,132L),2L),repeat_role=rep(c("minimum","median","maximum"),2L));primary<-mv17g_grouped_queue_v1("primary");repeat_queue<-mv17g_grouped_queue_v1("repeat",repeat_orders=repeat_orders);expect_equal(nrow(repeat_queue),27L);expect_equal(sum(repeat_queue$scientific_runs),2085L);expect_equal(sum(repeat_queue$null_family=="observed"),6L);keys<-c("view","unit_order","null_family","replicate_count","seed_first","seed_last");take<-paste(primary$view,primary$unit_order)%in%paste(repeat_orders$view,repeat_orders$unit_order);actual<-repeat_queue[keys];expected<-primary[take,keys];rownames(actual)<-rownames(expected)<-NULL;expect_equal(actual,expected)
})

test_that("MV17-G measured projection fits the prospective serial envelope", {
  resource<-data.frame(mode=c("primary","repeat"),jobs=c(195L,65L),outer_wall_seconds=c(816.83,296.40),private_bytes=c(5251136,4895357));p<-mv17g_resource_projection_v1(resource);expect_equal(p$grouped_children,c(1188L,27L));expect_equal(p$scientific_runs,c(91740L,2085L));expect_equal(unique(p$conservative_seconds_per_run),4.56);expect_equal(sum(p$projected_wall_seconds)/3600,118.845);expect_equal(sum(p$projected_wall_seconds_with_25pct_margin)/3600,148.55625);expect_true(all(p$projected_wall_seconds_with_25pct_margin<p$aggregate_timeout_seconds));expect_true(all(p$workers==1L&p$retries==0L))
})

test_that("MV17-G design freezes precision and every downstream firewall", {
  c<-mv17g_design_contract_v1();expect_equal(c$replicates,99L);expect_equal(c$minimum_attainable_probability,.01);expect_equal(c$worst_case_MCSE,.05);expect_true(c$H0_H1_separate);expect_false(c$full_calibration_authorized);expect_false(c$real_localization_authorized);expect_false(any(c[c("labels_authorized","outcomes_authorized","clustering_authorized","fusion_authorized","biology_authorized","manuscript_claims_authorized")]))
})
