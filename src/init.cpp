// #include <Rinternals.h>
// #include <R_ext/Rdynload.h>
// // #include "fitHMM_context.cpp"
// #include "RcppExports.cpp"
// 
// 
// R_NativePrimitiveArgType arg1[] = {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, VECSXP, INTSXP};
// 
// static const R_CMethodDef CEntries[]  = {
//     {"methimpute_fitBinomialTestHMMcontextTransition", (DL_FUNC) &fitBinomialTestHMMcontextTransition, 7, arg1},
// 		{"methimpute_cleanup", (DL_FUNC) &cleanup, 0, NULL},
//     {NULL, NULL, 0, NULL}
// };
// 
// 
// extern "C" {
// 	void R_init_methimpute(DllInfo *dll)
// 	{
// 		R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
// 		R_useDynamicSymbols(dll, FALSE);
// 	}
// }
