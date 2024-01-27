%============================================================================================
%    DFL - Derivative-Free Linesearch program for Mixed Integer Nonlinear Programming 
%    Copyright (C) 2011  G.Liuzzi, S.Lucidi, F.Rinaldi
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    G.Liuzzi, S.Lucidi, F.Rinaldi. Derivative-free methods for bound constrained 
%    mixed-integer optimization, Computational Optimization and Applications, 
%    53:505-526 (2012). DOI: 10.1007/s10589-011-9405-3
%============================================================================================
%
%-----------------------------------------------------------------------------!
%          Program DFL for bound-constrained mixed integer problems           !
%-----------------------------------------------------------------------------! 
%
%
%it solves the following problem:
%                                
%                             min f(x)
%                             l<=x<=u            
%                            
%                            l>-\infty
%                            u<\infty
%
%-----------------------------------------------------------------------------!
function [x,f,ni,nf,istop] = sd_box(funct,x0,index_int,scale_int,lb,ub,alfa_stop,nf_max,iprint)

    n  = length(x0);
    x  = x0;
    f  = feval(funct,x);
    nf = 1;
    ni = 0;
    
    discr_change = false; 
    eta          = 1.e-6;
    flag_fail    = false*ones(n,1);
    num_fal      = 0;
    istop        = 0;
    alfa_d       = zeros(n,1);
    
    %     ---- choice of the starting stepsizes along the directions ------
    alfa_max = 0.0;
    for i = 1:n
        if (index_int(i) == 0)
            alfa_d(i) = max(1.e-3,min(1.0,abs(x(i))));
            alfa_max  = max(alfa_max,alfa_d(i));
        else
            alfa_d(i) = max(scale_int(i),min(2.0*scale_int(i),abs(x(i))));
        end
        if (iprint >= 1)
            fprintf(' alfainiz(%3d) = %26.16e\n',i,alfa_d(i));
        end
    end
    
    d      = ones(n,1);
    i_corr = 1;
    z      = x;

    if (iprint >= 2)
        fprintf(' -------------------------------------------------\n');
        fprintf(' finiz = %26.16e\n',f);
        for i = 1:n
            fprintf(' xiniz(%3d) = %26.16e\n',i,x(i));
        end
        fprintf(' -------------------------------------------------\n');
    end

    % main loop of the algorithm
    while (true)
        if (iprint >= 0)
            fprintf(' ni=%4d  nf=%5d   f=%12.5e   alfamax=%12.5e\n',ni,nf,f,alfa_max); 
        end
        if (iprint >= 2)
            for i = 1:n
                fprintf(' x(%3d) = %26.16e\n',i,x(i));
            end
        end
        
%-------------------------------------
%    sampling along coordinate i_corr
%-------------------------------------
        if (index_int(i_corr) == 0)
            alfa = linesearchbox_cont();

        else
            alfa_d_old = alfa_d(i_corr);
            alfa = linesearchbox_discr();
        end
        
        if (abs(alfa) >= 1.e-12)
            flag_fail(i_corr) = false;
            x(i_corr)         = x(i_corr) + alfa*d(i_corr);
            f                 = fz;
            num_fal           = 0;
        else
            flag_fail(i_corr) = true;

            if ((index_int(i_corr) == 1) && (alfa_d_old > scale_int(i_corr)))  
                flag_fail(i_corr) = false;
            end
            if (i_corr_fall < 2)
                num_fal       = num_fal+1;
            end
        end

        ni        = ni+1;
        z(i_corr) = x(i_corr);

        if (i_corr < n)
            i_corr = i_corr+1;
        else
            if(~discr_change)
                for i = 1:n
                    if ((index_int(i) == 1) && (alfa_d(i) > 1))
                        discr_change = true;
                    	break
                    end
                end
                if (~discr_change)
                    eta = eta/2.0;
                end
            end
            i_corr       = 1;
            discr_change = false; 
        end

        istop = stop();

        if (istop >= 1)
            break
        end
    end
    if (iprint >= -1) 
        fprintf(' ni=%4d  nf=%5d   f=%12.5e   alfamax=%12.5e\n',ni,nf,f,alfa_max); 
    end

    function [istopout] = stop()
        alfa_max = 0.0;
        istopout = 0;
        for i = 1:n
            if (index_int(i) == 0)
                if (alfa_d(i) > alfa_max)
                    alfa_max = alfa_d(i);
                end
            end
        end
        
        if (alfa_max <= alfa_stop)
            test = true;
            for i = 1:n
                if (index_int(i) == 1)
                    if ((alfa_d(i) ~= scale_int(i)) || (flag_fail(i) ~= 1))
                        test = false;
                    end
                end
            end
            if (test)
                istopout = 1;
            end
        end
        
        if (nf > nf_max)
            istopout = 2;
        end

        if (ni > nf_max)
            istopout = 3;
        end
        
        return
    end

    function [alfa] = linesearchbox_cont()
        gamma  = 1.e-6;
        delta  = 0.5;
        delta1 = 0.5;
        ifront = 0;
        j      = i_corr;

        i_corr_fall = 0;

        if (iprint >= 1)
            fprintf('variabile continua  j =%3d    d(j) =%13.6e alfa=%13.6e\n',j,d(j),alfa_d(j));
        end
        
        if (abs(alfa_d(j)) <= 1.e-3*min(1.0,alfa_max))
            alfa = 0.0;
            if (iprint >= 1)
                fprintf('  alfa piccolo\n');
                fprintf(' alfa_d(j)=%13.6e    alfamax=%13.6e\n',alfa_d(j),alfa_max);
            end
            return
        end
        
        %choice of the direction
        for ielle = 1:2
            if (d(j) > 0.0)

                if ((alfa_d(j)-(ub(j)-x(j))) < -1.e-6)
                    alfa   = max(1.e-24,alfa_d(j));
                else
                    alfa   = ub(j)-x(j);
                    ifront = 1;
                    if (iprint >= 1)
                        fprintf(' point on the boundary. *\n');
                    end
                end

            else

                if ((alfa_d(j)-(x(j)-lb(j))) < -1.e-6)
                    alfa   = max(1.e-24,alfa_d(j));
                else
                    alfa   = x(j)-lb(j);
                    ifront = 1;
                    if (iprint >= 1)
                        fprintf(' point on the boundary. *\n');
                    end
                end

            end
            
            if (abs(alfa) <= 1.e-3*min(1.0,alfa_max))
                d(j)        = -d(j);
                i_corr_fall = i_corr_fall+1;
                ifront      = 0;

                if (iprint >= 1)
                    fprintf(' direzione opposta per alfa piccolo\n');
                    fprintf(' j =%d    d(j) =%e\n',j,d(j));
                    fprintf(' alfa=%e    alfamax=%e\n',alfa,alfa_max);
                end
                alfa = 0.0;
                continue
            end

            alfaex = alfa;
            z(j)   = x(j)+alfa*d(j);
            fz     = feval(funct,z);
            nf     = nf+1;
            fpar   = f - gamma*alfa*alfa;

            if (iprint >= 1)
                fprintf(' fz =%e   alfa =%e\n',fz,alfa);
            end
            if (iprint >= 2)
                for i = 1:n
                    fprintf(' z(%d)=%e\n',i,z(i));
                end
            end

            % test on the direction
            if (fz < fpar)
                % expansion step
                while (true)
                    if (ifront == 1)
                        if (iprint >= 1)
                            fprintf(' accetta punto sulla frontiera fz =%e   alfa =%e\n',fz,alfa);
                        end
                        alfa_d(j) = delta*alfa;
                        return
                    end

                    if (d(j) > 0.0)
                        if ((alfa/delta1-(ub(j)-x(j))) < -1.e-6)
                            alfaex = alfa/delta1;
                        else
                            alfaex = ub(j)-x(j);
                            ifront = 1;
                            if (iprint >= 1)
                                fprintf(' punto espan. sulla front.\n');
                            end
                        end

                    else
                        if ((alfa/delta1-(x(j)-lb(j))) < -1.e-6)
                            alfaex = alfa/delta1;
                        else
                            alfaex = x(j)-lb(j);
                            ifront = 1;
                            if (iprint >= 1)
                                fprintf(' punto espan. sulla front.\n');
                            end
                        end

                    end
                    z(j)    = x(j) + alfaex*d(j); 
                    fzdelta = feval(funct,z);
                    nf      = nf+1;
                    fpar    = f - gamma*alfaex*alfaex;

                    if (iprint >= 1)
                        fprintf(' fzex=%e  alfaex=%e\n',fzdelta,alfaex  );
                    end
                    if (iprint >= 2)
                        for i = 1:n
                            fprintf(' z(%d)=%e\n',i,z(i));
                        end
                    end

                    if (fzdelta < fpar)
                        fz   = fzdelta;
                        alfa = alfaex;
                    else               
                        alfa_d(j) = delta*alfa;
                        if (iprint >= 1)
                            fprintf(' accetta punto fz =%e   alfa =%e\n',fz,alfa);
                        end
                        return
                    end
                end % while(true)
            else
                d(j)   = -d(j);
                ifront = 0;
                if (iprint >= 1)
                    fprintf(' direzione opposta\n');
                    fprintf(' j =%d    d(j) =%e\n',j,d(j));
                end
            end % test on the direction
        end % for ielle = 1:2  
        
        if (i_corr_fall ~= 2)
            alfa_d(j) = delta*alfa_d(j);
        end
        
        alfa = 0.0;
        if (iprint >= 1)
            fprintf(' failure along the direction\n');
        end
        return      
    end

    function [alfa] = linesearchbox_discr()

        gamma_int   = 1.0; 
        delta       = 0.5;
        delta1      = 0.5;
        i_corr_fall = 0;
        ifront      = 0;
        j           = i_corr;
    
        if (iprint >= 1)
            fprintf('variabile discreta  j =%3d    d(j) =%13.6e alfa=%13.6e\n',j,d(j),alfa_d(j));
        end
        
        test = true;
        if (alfa_d(i_corr) == scale_int(i_corr))
            for i = 1:n
                if((flag_fail(i) == 0))
                    test = false;
                    break
                end
            end
            if (test)
                alfa = 0.0;
                if (iprint >= 1)
                    fprintf('direzione gia analizzata\n');
                end
                return
            end
        end
        
        for ielle = 1:2
            if (d(j) > 0.0)

                if ( ((ub(j)-x(j)-alfa_d(j) ) ) < 0.0 )
                    alfa   = ub(j)-x(j);
                    ifront = 1;
                    if (alfa == 0.0)
                        d(j)   = -d(j);
                        ifront = 0;
                        continue
                    end
                else
                    alfa = alfa_d(j);
                end

            else

                if ( ((x(j)-alfa_d(j)-lb(j))) < 0.0 )
                    alfa   = x(j)-lb(j);
                    ifront = 1;
                    if (alfa == 0.0)
                        d(j)   = -d(j);
                        ifront = 0;
                        continue
                    end
                else
                    alfa = alfa_d(j);
                end

            end

            alfaex = alfa;
            z(j)   = x(j)+alfa*d(j);
            fz     = feval(funct,z);
            nf     = nf+1;
            fpar   = f-gamma_int*eta;
                        
            if (iprint >= 1)
                fprintf(' fz =%e   alfa =%e\n',fz,alfa);
            end
            if (iprint >= 2)
                for i = 1:n
                    fprintf(' z(%d)=%e\n',i,z(i));
                end
            end

            if (fz < fpar)
                discr_change = true;
                while (true)
                    if (ifront == 1)
                        if (iprint >= 1)
                            fprintf(' accetta punto sulla frontiera fz =%e   alfa =%e\n',fz,alfa);
                        end
                        return
                    end

                    if (d(j) > 0.0)
                        if ((ub(j)-x(j)-2.0*alfa ) < 0.0)
                            alfaex = ub(j)-x(j);
                            ifront = 1;
                            if (alfaex <= alfa)
                                alfa_d(j) = max(scale_int(j),max(floor((alfa/2.0)/scale_int(j)+0.5),1.0)*scale_int(j));
                                if (iprint >= 1)
                                    fprintf(' accetta punto quasi frontiera fz =%e   alfa =%e\n',fz,alfa);
                                end
                                return
                            end
                        else
                            alfaex = alfa*2.0;
                        end
                    else
                        if (( x(j)-2.0*alfa-lb(j) ) < 0.0)
                            alfaex = x(j)-lb(j);
                            ifront = 1;
                            if (alfaex <= alfa)
                                alfa_d(j) = max(scale_int(j),max(floor((alfa/2.0)/scale_int(j)+0.5),1.0)*scale_int(j));
                                if (iprint >= 1)
                                    fprintf(' accetta punto quasi frontiera fz =%e   alfa =%e\n',fz,alfa);
                                end
                                return
                            end
                        else
                          alfaex = alfa*2.0;
                        end
                    end
                    
                    z(j)    = x(j) + alfaex*d(j); 
                    fzdelta = feval(funct,z);
                    nf      = nf+1;
                    fpar    = f - gamma_int*eta;

                    if (iprint >= 1)
                        fprintf(' fzex=%e  alfaex=%e\n',fzdelta,alfaex  );
                    end
                    if (iprint >= 2)
                        for i = 1:n
                            fprintf(' z(%d)=%e\n',i,z(i));
                        end
                    end

                    if (fzdelta < fpar)
                        fz   = fzdelta;
                        alfa = alfaex;
                    else               
                        alfa_d(j) = max(scale_int(j),max(floor((alfa/2.0)/scale_int(j)+0.5),1.0)*scale_int(j));
                        if (iprint >= 1)
                            fprintf(' accetta punto fz =%e   alfa =%e\n',fz,alfa);
                        end
                        return
                    end
                end % while (true)
            else
                d(j)   = -d(j);
                ifront = 0;

                if (iprint >= 1)
                    fprintf(' direzione opposta\n');
                    fprintf(' j =%d    d(j) =%e\n',j,d(j));
                end

            end % if (fz < fpar)
        end
        alfa_d(j) = max(scale_int(j),max(floor((alfa/2.0)/scale_int(j)+0.5),1.0)*scale_int(j));
        alfa      = 0.0;
        if (iprint >= 1)
            fprintf(' failure along the direction\n');
        end
        return
    end
end
