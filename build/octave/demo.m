my_benchmark = benchmark("bbob", "bbob", "random_search")
while (true)
    problem = next_problem(my_benchmark)
    my_optimizer(problem, problem.lower_bounds, problem.upper_bounds, 1000)
    free(problem)
endwhile