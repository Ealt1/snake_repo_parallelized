classdef SASys<handle

	properties (SetAccess=protected)
        A;
        K;
        t; 
	end

	methods
		function sys = SASys(A,K,t)
			% SASys: Create a instant of Fuel Cell Thermal Managment
			%          System, under one specific control u
			% 
			% SYNTAX
			% ------
			%
			%	
			%
			% INPUT
			% -----
			%	
			%	
			%
			% PROPERTIES
			% ------
			%
            %   see INPUT
			% 
			% SEE ALSO
			% -------
			%
            sys.A = A;
            sys.K = K;
            sys.t = t;
        end
        
        function r = evalVField(sys, xx, idx)
            % evalVField: 
            % Input: 
            %   sys: SASys
            %   x: [x1,...,xm]
            %   idx: index set, subset of [1:7]
            % Output: 
            %   r: [x+_{idx}1,...,x+_{idx}m]
            if (nargin < 3)
                n = size(sys.A,1);
                idx = 1:n;
            end
            r = [];
            for i = 1:1:size(xx,2)
                x = xx(:,i);
                vfields = sys.A*x + sys.K;
                r = [r, vfields(idx)];
            end
        end
        
        function X_plus = evolve_Rec(sys,X)
            % input: 
            %   X, Rec
            % output: 
            %   X_plus, is the smallest Rec containing set 
            %   {sys.A*x + sys.K: x \in X}
            lux = [X.xmin', X.xmax'];
            m = [];
            M = [];
            for i = 1:1:X.dim
                sppx = sys.A(i,:)';
                % minimal (or minimizer)
                mx_idx = (sppx<=0)+1;
                Mx_idx = 3 - mx_idx;

                mx = diag(lux(:,mx_idx));
                Mx = diag(lux(:,Mx_idx));

                m_i = min(evalVField(sys, mx, i));
                M_i = max(evalVField(sys, Mx, i));
                
                m = [m, m_i];
                M = [M, M_i];
            end
            X_plus = Rec([m; M]);
        end
        
        function X_plus = evolve_Poly(sys,X)
            % input: X, Polyhedron
            % output: X_plus, Polyhedron
            V_plus = sys.A*X.V' + repmat(sys.K,1,size(X.V,1));
            X_plus = Polyhedron('V',V_plus');
        end
        
        function X_plus = evolve_Rec_by_Poly(sys,X)
            % input: X, Polyhedron
            % output: X_plus, Polyhedron
            P = Polyfy(X); 
            V_plus = sys.A*P.V' + repmat(sys.K,1,size(P.V,1));
            P_plus = Polyhedron('V',V_plus');
            X_plus = Rec([min(P_plus.V,[],1); max(P_plus.V,[],1)]);
        end
	end
end


