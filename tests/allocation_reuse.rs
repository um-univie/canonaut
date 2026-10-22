//! Count requests on the test thread, excluding unrelated test-runner traffic.
use canonaut::{generators::generate_random_graph, structs::{CanonautManager, ValidatedGraph}};
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

thread_local! {
    static TRACK: Cell<bool> = const { Cell::new(false) };
    static REQUESTS: Cell<usize> = const { Cell::new(0) };
}
struct CountingAllocator;
fn count() {
    let _ = TRACK.try_with(|track| {
        if track.get() {
            let _ = REQUESTS.try_with(|n| n.set(n.get() + 1));
        }
    });
}
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        count();
        unsafe { System.alloc(layout) }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        count();
        unsafe { System.alloc_zeroed(layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        count();
        unsafe { System.realloc(ptr, layout, size) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
}
#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

#[test]
fn warmed_colored_and_uncolored_size_changes_do_not_allocate() {
    let graph = generate_random_graph(193, 0.1, 42);
    let small = vec![0; 127 * 2];
    let colors: Vec<u32> = (0..127).map(|v| if v < 3 { 0 } else { v }).collect();
    let mut manager = CanonautManager::new(193).with_canonization();
    manager.canonize(&graph);
    manager.canonize_with_colors(&small, &colors);
    let validated = ValidatedGraph::new(&graph);
    let validated_small = ValidatedGraph::new(&small);
    REQUESTS.with(|n| n.set(0));
    TRACK.with(|track| track.set(true));
    for _ in 0..20 {
        manager.canonize(&graph);
        manager.canonize_with_colors(&small, &colors);
        manager.canonize_validated(&validated);
        manager.canonize_validated_with_colors(&validated_small, &colors);
    }
    TRACK.with(|track| track.set(false));
    assert_eq!(REQUESTS.with(Cell::get), 0);
    assert_eq!(manager.statistics().group_size1, 6.0);
}
