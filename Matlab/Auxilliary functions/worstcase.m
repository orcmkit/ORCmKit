function [x_LB,x_UB,f_LB,f_MID,f_UB,minus_percent,plus_percent] = worstcase(f,x,uncertainties)
%[X_LB,X_UB,F_LB,F_MID,F_UB] = WORSTCASE(F,X,UNCERTAINTIES)
%
%Performs worst-case maximization of the error realized in a function f,
%about a vector x, when there is an error associated with each element in x.
%worstcase is a different (and much better) method of calculating the
%propagation of uncertainty in a function. This method is superior to the
%prior function I wrote ("Propagation of Uncertainty", or PropError.m)for
%several reasons:
%
%i.) it's faster,
%ii.) greater input flexibility
%iii.) does not require symbolic functions or the symbolic math toolbox; it
%only requires the more-commonplace optimization toolbox,
%iv.) does not require the user to put in a list of the variables
%v.) works with implicit functions, and functions which lack an analytic
%derivative.
%
%worstcase solves the following optimization problems using fmincon:
%min f(d), subject to x(i) - uncertainties(i) <= d(i) <= x(i) + uncertainties(i).
%max f(d), subject to x(i) - uncertainties(i) <= d(i) <= x(i) + uncertainties(i).
%
%The default initial guess x0 is at the vector x.
%
%f is a function that accepts a vector, x, and outputs a scalar value. If x
%is a matrix, then worstcase finds the worst-case for each column in x.
%
%uncertainties is a vector of absolute uncertainties in the values of x. All 
%elements of uncertainties must be positive. uncertainties can be either a single column
%or a matrix. If uncertainties is a single column, that error vector is applied to
%all of the columns in x. If uncertainties is a matrix, the columns of uncertainties are
%considered to be the uncertainties in a corresponding column in x. e.g.,
%the optimization problem will be solved for x(:,i) and uncertainties(:,i).
%
%The default algorithm is sqp. If f is highly nonconvex, or if uncertainties(i) is
%large relative to x(i), sqp may become trapped in a suboptimal local
%minimum. If the exitflag for sqp is not 1 and not 2, then a warning is
%displayed to alert the user to sub-optimality.
%
%The outputs of worstcase are:
%x_LB, the "lower bound" x-matrix. The columns of x_LB produced the
%minimized f-values in f_LB. f_LB is the lower bound for the function, f.
%x_UB and f_UB are corresponding matrices for the upper bound. f_MID is
%simply f(x); the value of f for each column of x with zero error.
%minus_percent and plus_percent are conversions of the calculated
%uncertainty in f to percentage form.
%
%(c) Bradley James Ridder, October 15, 2014. Purdue University, West
%Lafayette, Indiana, USA.

    myopts = optimset('fmincon');
    myopts = optimset(myopts,'Display','off','Algorithm','sqp');

    nvec = size(x,2);
    nvars = size(x,1);

    if size(uncertainties,1) ~= nvars
        error('Number of rows in x must equal the number of rows in uncertainties.')
    end
    
    if any(uncertainties(:) < 0)
        error('uncertainties must all be positive.')
    end
    
    x_LB = zeros(nvars,nvec);
    x_UB = zeros(nvars,nvec);
    f_LB = zeros(    1,nvec);
    f_UB = zeros(    1,nvec);
    f_MID = zeros(   1,nvec);
    exitflag_LB = zeros(1,nvec);
    exitflag_UB = zeros(1,nvec);
    for i = 1:nvec
        if size(uncertainties,2) > 1
            if size(uncertainties,2) ~= nvec
                error('If uncertainties has more than one column, the number of columns must equal the number of columns in x.')
            end
            q = uncertainties(:,i);
        else
            q = uncertainties;
        end
        f_MID(i) = f(x(:,i));
        [x_LB(:,i),f_LB(i),exitflag_LB(1,i)] = fmincon(@(x)         f(x) ,x(:,i),[],[],[],[],x(:,i) - q,x(:,i) + q,[],myopts);
        [x_UB(:,i),f_UB(i),exitflag_UB(1,i)] = fmincon(@(x) minus(0,f(x)),x(:,i),[],[],[],[],x(:,i) - q,x(:,i) + q,[],myopts);
        if exitflag_LB(1,i) ~= 1 && exitflag_LB(1,i) ~= 2
            warning('Possible sub-optimal solution found when calculating worst-case LOWER bound; x-column %d, exitflag %d.',i,exitflag_LB(1,i))
        end
        if exitflag_UB(1,i) ~= 1 && exitflag_UB(1,i) ~= 2
            warning('Possible sub-optimal solution found when calculating worst-case UPPER bound; x-column %d, exitflag %d.',i,exitflag_UB(1,i))
        end
    end
    f_UB = -f_UB; %We took the negative before when we were maximizing it.
                   %Now we will take the negative again to put ourselves back on our feet.
    minus_percent = strcat(num2str(-abs((1- f_LB./f_MID)*100)','%0.5g'),'%');
    plus_percent =  strcat('+',num2str( abs((1- f_UB./f_MID)*100)','%0.5g'),'%');    
end