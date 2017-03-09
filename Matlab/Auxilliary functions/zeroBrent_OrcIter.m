function value = zeroBrent_OrcIter ( i, a, b, machep, t, ti, f )

%*****************************************************************************80
%
%% ZERO seeks the root of a function F(X) in an interval [A,B].
%
%  Discussion:
%
%    The interval [A,B] must be a change of sign interval for F.
%    That is, F(A) and F(B) must be of opposite signs.  Then
%    assuming that F is continuous implies the existence of at least
%    one value C between A and B for which F(C) = 0.
%
%    The location of the zero is determined to within an accuracy
%    of 6 * MACHEPS * abs ( C ) + 2 * T.
%
%    Thanks to Thomas Secretin for pointing out a transcription error in the
%    setting of the value of P, 11 February 2013.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    11 February 2013
%
%  Author:
%
%    Original FORTRAN77 version by Richard Brent
%    MATLAB version by John Burkardt
%
%  Reference:
%
%    Richard Brent,
%    Algorithms for Minimization Without Derivatives,
%    Dover, 2002,
%    ISBN: 0-486-41998-3,
%    LC: QA402.5.B74.
%
%  Parameters:
%
%    Input, real A, B, the endpoints of the change of sign interval.
%
%    Input, real MACHEP, an estimate for the relative machine
%    precision.
%
%    ti, residuals on the guess residuals
%
%    Input, real T, a positive error tolerance.
%
%    Input, real value = F ( x ), the name of a user-supplied
%    function which evaluates the function whose zero is being sought.
%
%    Output, real VALUE, the estimated value of a zero of
%    the function F.
%

%
%  Make local copies of A and B.

% Custom Rémi DICKES 08/02/2017
% If there is a guess on the subcooling provided by the user
if not(isnan(i))
    fi = f(i);
    if abs(fi)<ti
        value = i;
        return
    end
    
    if i>0 
        delta = 5;
        if fi >0
            i_delta = min(i + delta, b);
            fi_delta = f(i_delta);
        else
            i_delta = max(i - delta, 0);
            fi_delta = f(i_delta);            
        end

    else
        delta = 1e-3;
        if fi >0
            i_delta = min(i + delta, 0);
            fi_delta = f(i_delta);
        else
            i_delta = max(i - delta, a);
            fi_delta = f(i_delta);            
        end

    end
    
    if abs(fi_delta)<ti
        value = i_delta;
        return
    end
    
else
    fi = NaN;
    i_delta = NaN;
    fi_delta = NaN;
end


if fi > 0 && fi_delta < 0
    sa = i;
    fa = fi;
    sb = i_delta;
    fb = fi_delta;
elseif fi < 0 && fi_delta > 0
    sa = i_delta;
    fa = fi_delta;
    sb = i;
    fb = fi;    
    
else
    i = i_delta;
    fi = fi_delta;
        
    f0 = f(0);
    if f0 == 0
        value = 0;
        return
        
    elseif f0 > 0
        sa = 0;
        fa = f0;
        if fi < 0
            sb = i;
            fb = fi;
        else
            sb = b;
            fb = f ( sb );
            if fb == 0
                value = sb;
                return
            elseif fb > 0
                sb = 2*b;
                fb = f(sb);
                if fb > 0 
                    value = sb;
                    return
                end
            end
        end
        
    elseif f0 < 0
        sb = 0;
        fb = f0;
        if fi > 0
            sa = i;
            fa = fi;
        else
            sa = a;
            fa = f(sa);
            if fa == 0
                value = sa;
                return
            elseif fa < 0
                sa = 2*a;
                fa = f(sa);
                if fa < 0
                    value = sa;
                    return
                end
            end
        end
    end
end

c = sa;
fc = fa;
e = sb - sa;
d = e;

while ( 1 )
    
    if ( abs ( fc ) < abs ( fb ) )
        
        sa = sb;
        sb = c;
        c = sa;
        fa = fb;
        fb = fc;
        fc = fa;
        
    end
    
    tol = 2.0 * machep * abs ( sb ) + t;
    m = 0.5 * ( c - sb );
    
    if ( abs ( m ) <= tol || fb == 0.0 || abs(fb) <= ti)
        break
    end
    
    if ( abs ( e ) < tol || abs ( fa ) <= abs ( fb ) )
        
        e = m;
        d = e;
        
    else
        
        s = fb / fa;
        
        if ( sa == c )
            
            p = 2.0 * m * s;
            q = 1.0 - s;
            
        else
            
            q = fa / fc;
            r = fb / fc;
            p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
            q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
            
        end
        
        if ( 0.0 < p )
            q = - q;
        else
            p = - p;
        end
        
        s = e;
        e = d;
        
        if ( 2.0 * p < 3.0 * m * q - abs ( tol * q ) && p < abs ( 0.5 * s * q ) )
            d = p / q;
        else
            e = m;
            d = e;
        end
        
    end
    
    sa = sb;
    fa = fb;
    
    if ( tol < abs ( d ) )
        sb = sb + d;
    elseif ( 0.0 < m )
        sb = sb + tol;
    else
        sb = sb - tol;
    end
    
    fb = f ( sb );
    
    if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
        c = sa;
        fc = fa;
        e = sb - sa;
        d = e;
    end
    
end

value = sb;

return
end