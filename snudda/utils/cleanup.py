def cleanup_memory(d_view, variable_list):
    """ Clear variables in variable_list on workers specified by d_view.

    Args:
        d_view : ipyparallel direct view
        variable_list : list of variables to clear
        """

    for var in variable_list:
        d_view.execute(f"if '{var}' in locals(): del {var}")
        # d_view.execute(f"{var}=None")


# Cleans up memory on workers after execution
def cleanup(rc, state):
    """ When place, detect, prune are run one after the other, we do not need to keep them all in memory
    at the same time. As one finishes, we can clear the memory.

    Args:
        rc : ipyparallel remote client
        state : State to free memory (e.g. "place", "detect", "project", "prune"
        """

    if rc is not None:

        # Check if the client is still active
        try:
            if hasattr(rc, '_closed') and rc._closed:
                return
            # Alternative check: try to access a property that would fail if closed
            _ = rc.ids
        except (RuntimeError, AttributeError):
            # Client is closed or invalid, skip cleanup
            return

        d_view = rc.direct_view(targets='all')

        var_lookup = {"place": ["inner_mask", "sm"],  # region_mesh_vedo.py
                      "detect": ["min_max", "result", "sd"],
                      "project": [],  # Currently does not support parallel execution
                      "prune": ["syn_before", "syn_after", "merge_result_syn", "merge_result_gj", "sp"]}

        if state in var_lookup:
            cleanup_memory(d_view, var_lookup[state])
        else:
            print(f"cleanup: Unknown state: {state}")
