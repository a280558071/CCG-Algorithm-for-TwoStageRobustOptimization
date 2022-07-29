# CCG Algorithm for Two-Stage Robust Optimization
 This project contains two MATLAB files explain Benders-dual/Column and Constraints Generation (CCG) Algorithm to solve Two-Stage Robust Optimization (RO) problems. CCG is proposed in [1]. It's very fashion in optimization for power systems/integrated energy systems and has been cited 1000+ times.
 - BDforExample.m explains benders-dual method to solve Two-Stage RO problems
 - CCGforExample.m explains CCG Algorithm to solve Two-Stage RO problems
 - The example used in both files are from *4. Case study: robust location-transportation problem* in [1]
 - MATLAB optimization modeling tool YALMIP supports constraints writing in the fashion of “Cons = [A*x<=b]”, in which A is a matrix and x, b are vectors, therefore, the codes in both files are easy to read. (at least I believe so...)
 - Because they are just like what you can see about the mathematical model in the paper.
# Must-include package
 - For both files: **YALMIP**, see (https://yalmip.github.io/)
 - Both files are calling GUROBI(https://www.gurobi.com/) to solve the problem, but other solvers (e.g. CPLEX) could be used as well. (see what solvers they support in abovementioned link for YALMIP).
# What you (maybe) can learn
 - Basics of Benders-dual/CCG Algorithm to solve Two-Stage RO problems
 - A simple way to realize optimization modeling in MATLAB (ONLY ~100 lines are in both .m files)
# What you CANNOT learn
 - How to do in by yourself (just try it by yourself and you would definitely learn more!)
# Acknowledgement
 - [鲁棒优化| C&CG算法求解两阶段鲁棒优化：全网最完整、最详细的【入门-完整推导-代码实现】笔记](https://zhuanlan.zhihu.com/p/534285185). I found some tricky point in it about how to guarantee Master problem could provide initial feasible solution for Subproblem to make both algorithms work, so thanks.
# To be updated/some issue
 - When Subproblem is unbounded both file would fail. Because I still don't know how to identify scenario for which Q(y*)=inf.
# Ref
[1] Zeng, Bo, and Long Zhao. "Solving two-stage robust optimization problems using a column-and-constraint generation method." Operations Research Letters 41, no. 5 (2013): 457-461.
