def cleanup_memory(d_view, variable_list):

    for var in variable_list:
        d_view.execute(f"if '{var}' in locals(): del {var}")
        # d_view.execute(f"{var}=None")


# Cleans up memory on workers after execution
def cleanup(rc, state):

    d_view = rc.direct_view(targets='all')

    var_lookup = {"place": ["innerMask", "sm"],  # region_mesh.py
                  "detect": ["min_max", "result", "sd"],
                  "project": [],  # Currently does not support parallel execution
                  "prune": ["syn_before", "syn_after", "merge_result_syn", "merge_result_gj", "sp"]}

    if state in var_lookup:
        cleanup_memory(d_view, var_lookup[state])
    else:
        print(f"cleanup: Unknown state: {state}")
