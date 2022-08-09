#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/*
 * Required to make our bridge to demes-forward-capi 
 * know that this type will be provided at link-time
 */
typedef struct OpaqueForwardGraph OpaqueForwardGraph;

#include "internal/fp11rust_api.h"

#ifdef __cplusplus
} // extern "C"
#endif // __cplusplus
