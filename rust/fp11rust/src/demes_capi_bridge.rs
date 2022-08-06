use demes_forward_capi::*;
use libc::c_char;

#[no_mangle]
pub extern "C" fn demes_forward_graph_allocate() -> *mut OpaqueForwardGraph {
    forward_graph_allocate()
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_initialize_from_yaml(
    yaml: *const c_char,
    burnin: f64,
    graph: *mut OpaqueForwardGraph,
) -> i32 {
    forward_graph_initialize_from_yaml(yaml, burnin, graph)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_initialize_from_yaml_file(
    file_name: *const c_char,
    burnin: f64,
    graph: *mut OpaqueForwardGraph,
) -> i32 {
    forward_graph_initialize_from_yaml_file(file_name, burnin, graph)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_is_error_state(
    graph: *const OpaqueForwardGraph,
) -> bool {
    forward_graph_is_error_state(graph)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_deallocate(graph: *mut OpaqueForwardGraph) {
    forward_graph_deallocate(graph)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_get_error_message(
    graph: *const OpaqueForwardGraph,
    status: *mut i32,
) -> *const c_char {
    forward_graph_get_error_message(graph, status)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_selfing_rates(
    graph: *const OpaqueForwardGraph,
    status: *mut i32,
) -> *const f64 {
    forward_graph_selfing_rates(graph, status)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_cloning_rates(
    graph: *const OpaqueForwardGraph,
    status: *mut i32,
) -> *const f64 {
    forward_graph_cloning_rates(graph, status)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_parental_deme_sizes(
    graph: *const OpaqueForwardGraph,
    status: *mut i32,
) -> *const f64 {
    forward_graph_parental_deme_sizes(graph, status)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_offspring_deme_sizes(
    graph: *const OpaqueForwardGraph,
    status: *mut i32,
) -> *const f64 {
    forward_graph_offspring_deme_sizes(graph, status)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_any_extant_offspring_demes(
    graph: *const OpaqueForwardGraph,
    status: *mut i32,
) -> bool {
    forward_graph_any_extant_offspring_demes(graph, status)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_any_extant_parent_demes(
    graph: *const OpaqueForwardGraph,
    status: *mut i32,
) -> bool {
    forward_graph_any_extant_parent_demes(graph, status)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_number_of_demes(
    graph: *const OpaqueForwardGraph,
) -> isize {
    forward_graph_number_of_demes(graph)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_update_state(
    time: f64,
    graph: *mut OpaqueForwardGraph,
) -> i32 {
    forward_graph_update_state(time, graph)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_initialize_time_iteration(
    graph: *mut OpaqueForwardGraph,
) -> i32 {
    forward_graph_initialize_time_iteration(graph)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_iterate_time(
    graph: *mut OpaqueForwardGraph,
    status: *mut i32,
) -> *const f64 {
    forward_graph_iterate_time(graph, status)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_ancestry_proportions(
    offspring_deme: usize,
    status: *mut i32,
    graph: *mut OpaqueForwardGraph,
) -> *const f64 {
    forward_graph_ancestry_proportions(offspring_deme, status, graph)
}

#[no_mangle]
pub unsafe extern "C" fn demes_forward_graph_model_end_time(
    status: *mut i32,
    graph: *const OpaqueForwardGraph,
) -> f64 {
    forward_graph_model_end_time(status, graph)
}
