import precapitate
import recapitate
%load_ext line_profiler
%lprun -f recapitate.run_sim recapitate.run_sim(1000, 1000, 1000, 1000, "dtwf", 100, 101, 666)
%lprun -f precapitate.run_sim precapitate.run_sim(1000, 1000, 1000, 1000, "dtwf", 100, 101, 666)
