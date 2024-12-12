; ModuleID = 'interSwap.bc'
source_filename = "../../../../TestCases/test/interSwap.c"
target datalayout = "e-m:o-i64:64-i128:128-n32:64-S128"
target triple = "arm64-apple-macosx15.0.0"

; Function Attrs: mustprogress nofree nosync nounwind ssp uwtable willreturn writeonly
define i32 @main() local_unnamed_addr #0 {
entry:
  %A.sroa.0 = alloca i8, align 1
  %B = alloca i8, align 1
  call void @llvm.lifetime.start.p0i8(i64 1, i8* nonnull %A.sroa.0)
  call void @llvm.lifetime.start.p0i8(i64 1, i8* nonnull %B)
  %A.sroa.0.0.A.sroa.0.0.A.sroa.0.0. = load i8, i8* %A.sroa.0, align 1, !tbaa !10
  %tobool.not = icmp eq i8 %A.sroa.0.0.A.sroa.0.0.A.sroa.0.0., 0
  br i1 %tobool.not, label %if.else, label %if.then

if.then:                                          ; preds = %entry
  store i8 65, i8* %A.sroa.0, align 1, !tbaa !10
  br label %if.end

if.else:                                          ; preds = %entry
  store i8 66, i8* %B, align 1, !tbaa !10
  br label %if.end

if.end:                                           ; preds = %if.else, %if.then
  %a.0 = phi i8* [ %B, %if.else ], [ %A.sroa.0, %if.then ]
  store i8 63, i8* %a.0, align 1, !tbaa !10
  call void @llvm.lifetime.end.p0i8(i64 1, i8* nonnull %B)
  call void @llvm.lifetime.end.p0i8(i64 1, i8* nonnull %A.sroa.0)
  ret i32 0
}

; Function Attrs: argmemonly nofree nosync nounwind willreturn
declare void @llvm.lifetime.start.p0i8(i64 immarg, i8* nocapture) #1

; Function Attrs: mustprogress nofree norecurse nosync nounwind ssp uwtable willreturn
define void @swap(i8** nocapture noundef %p, i8** nocapture noundef %q) local_unnamed_addr #2 {
entry:
  %0 = load i8*, i8** %p, align 8, !tbaa !13
  %1 = load i8*, i8** %q, align 8, !tbaa !13
  store i8* %1, i8** %p, align 8, !tbaa !13
  store i8* %0, i8** %q, align 8, !tbaa !13
  ret void
}

; Function Attrs: argmemonly nofree nosync nounwind willreturn
declare void @llvm.lifetime.end.p0i8(i64 immarg, i8* nocapture) #1

attributes #0 = { mustprogress nofree nosync nounwind ssp uwtable willreturn writeonly "frame-pointer"="non-leaf" "min-legal-vector-width"="0" "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="apple-m1" "target-features"="+aes,+crc,+crypto,+dotprod,+fp-armv8,+fp16fml,+fullfp16,+lse,+neon,+ras,+rcpc,+rdm,+sha2,+v8.5a,+zcm,+zcz" }
attributes #1 = { argmemonly nofree nosync nounwind willreturn }
attributes #2 = { mustprogress nofree norecurse nosync nounwind ssp uwtable willreturn "frame-pointer"="non-leaf" "min-legal-vector-width"="0" "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="apple-m1" "target-features"="+aes,+crc,+crypto,+dotprod,+fp-armv8,+fp16fml,+fullfp16,+lse,+neon,+ras,+rcpc,+rdm,+sha2,+v8.5a,+zcm,+zcz" }

!llvm.module.flags = !{!0, !1, !2, !3, !4, !5, !6, !7, !8}
!llvm.ident = !{!9}

!0 = !{i32 2, !"SDK Version", [2 x i32] [i32 15, i32 1]}
!1 = !{i32 1, !"wchar_size", i32 4}
!2 = !{i32 1, !"branch-target-enforcement", i32 0}
!3 = !{i32 1, !"sign-return-address", i32 0}
!4 = !{i32 1, !"sign-return-address-all", i32 0}
!5 = !{i32 1, !"sign-return-address-with-bkey", i32 0}
!6 = !{i32 7, !"PIC Level", i32 2}
!7 = !{i32 7, !"uwtable", i32 1}
!8 = !{i32 7, !"frame-pointer", i32 1}
!9 = !{!"clang version 14.0.6 (git@github.com:themaplelab/llvm-pointer-analysis.git 09cb671f4fc7974a29608cf9f122cf714461c4b2)"}
!10 = !{!11, !11, i64 0}
!11 = !{!"omnipotent char", !12, i64 0}
!12 = !{!"Simple C/C++ TBAA"}
!13 = !{!14, !14, i64 0}
!14 = !{!"any pointer", !11, i64 0}
