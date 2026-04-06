function objective = make_counted_objective(fun)
    n_calls = 0;

    function y = eval_fun(x)
        n_calls = n_calls + 1;
        y = fun(x);
    end

    function n = get_count_fun()
        n = n_calls;
    end

    objective = struct();
    objective.eval = @eval_fun;
    objective.get_count = @get_count_fun;
end