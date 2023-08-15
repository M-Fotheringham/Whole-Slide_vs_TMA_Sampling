###
def eta_counter(start_time, total_sampled, total_iterations):
    time_per_iteration = (datetime.datetime.now() - start_time) / len(total_sampled)
    time_to_go = time_per_iteration * (total_iterations - len(total_sampled))
    eta = (datetime.datetime.now() + time_to_go).strftime("%B %d, %Y %H:%M:%S")
    return eta, time_per_iteration
