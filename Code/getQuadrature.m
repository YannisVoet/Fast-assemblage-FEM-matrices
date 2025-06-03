function [points, weights] = getQuadrature(id, object)

% getQuadrature: Function which returns the integration points and
% weights for a given quadrature ID
% INPUT:
% id:       Integration ID
% OUTPUT:
% points:   Integration points
% weights:  Integration weights

switch object
    case 'bulk'
        switch id
            case 1
                % Source: Wriggers, Computational contact mechanics (2006)
                % Number of points: 1
                % Degree of exactness: 1
                
                a1=1/3;                 w1=1;
                
                points=[a1 a1];
                
                weights=0.5*w1;
                
            case 2
                % Source: Wriggers, Computational contact mechanics (2006)
                % Number of points: 3
                % Degree of exactness: 2
                
                a1=2/3;                 w1=1/3;
                b1=1/6;
                
                points=[b1 b1;
                        a1 b1;
                        b1 a1];
                
                weights=0.5*[w1 w1 w1];
                
            case 3
                % Source: Wriggers, Computational contact mechanics (2006)
                % Number of points: 4
                % Degree of exactness: 3
                
                a1=1/3;                 w1=-27/48;
                a2=3/5;                 w2=25/48;
                b2=1/5;
                
                points=[a1 a1;            
                        b2 b2;
                        a2 b2;
                        b2 a2];
                
                weights=0.5*[w1 w2 w2 w2];
                
            case 4
                % Source: Laursen and Gellert, Some criteria for
                % numerically integrated matrices and quadrature formulas
                % for triangles (1978)
                % Number of points: 6
                % Degree of exactness: 4
                
                a1=0.816847572980459;       w1=0.109951743655322;
                b1=0.091576213509771;       w2=0.223381589678011;
                a2=0.108103018168070;
                b2=0.445948490915965;
                
                points=[b1 b1;
                        a1 b1;
                        b1 a1;
                        b2 b2;
                        a2 b2;
                        b2 a2];
                
                weights=0.5*[w1 w1 w1 w2 w2 w2];
                
            case 5
                % Source: Laursen and Gellert, Some criteria for
                % numerically integrated matrices and quadrature formulas
                % for triangles (1978)
                % Number of points: 7
                % Degree of exactness: 5
                
                a1=1/3;                     w1=0.225000000000000;
                a2=0.797426985353087;       w2=0.125939180544827;
                b2=0.101286507323456;       w3=0.132394152788506;
                a3=0.059715871789770;
                b3=0.470142064105115;
                
                points=[a1 a1;
                        b2 b2;
                        a2 b2;
                        b2 a2;
                        b3 b3;
                        a3 b3;
                        b3 a3];
                    
               weights=0.5*[w1 w2 w2 w2 w3 w3 w3];
                    
            case 6
                % Source: Laursen and Gellert, Some criteria for
                % numerically integrated matrices and quadrature formulas
                % for triangles (1978)
                % Number of points: 12
                % Degree of exactness: 6
                
                a1=0.873821971016996;       w1=0.050844906370207;
                b1=0.063089014491502;       w2=0.116786275726379;
                a2=0.501426509658179;       w3=0.082851075618374;
                b2=0.249286745170910;
                a3=0.636502499121399;
                b3=0.310352451033785;
                c3=0.053145049844816;
                
                points=[b1 b1;
                        b1 a1;
                        a1 b1;
                        b2 b2;
                        b2 a2;
                        a2 b2;
                        a3 b3;
                        a3 c3;
                        b3 a3;
                        b3 c3;
                        c3 a3;
                        c3 b3];
                    
               weights=0.5*[w1 w1 w1 w2 w2 w2 w3 w3 w3 w3 w3 w3];
                 
            case 7
                % Source: Laursen and Gellert, Some criteria for
                % numerically integrated matrices and quadrature formulas
                % for triangles (1978)
                % Number of points: 13
                % Degree of exactness: 7
                
                a1=1/3;                     w1=-0.149570044467670;
                a2=0.479308067841923;       w2= 0.175615257433204;           
                b2=0.260345966079038;       w3= 0.053347235608839;
                a3=0.869739794195568;       w4= 0.077113760890257;
                b3=0.065130102902216;
                a4=0.638444188569809;
                b4=0.312865496004875;
                c4=0.048690315425316;
                
                points=[a1 a1;
                        b2 b2;
                        b2 a2;
                        a2 b2;
                        b3 b3;
                        b3 a3;
                        a3 b3;
                        a4 b4;
                        a4 c4;
                        b4 a4;
                        b4 c4;
                        c4 a4;
                        c4 b4];
                                                      
                weights=0.5*[w1 w2 w2 w2 w3 w3 w3 w4 w4 w4 w4 w4 w4];

            case 8
                % Source: Laursen and Gellert, Some criteria for
                % numerically integrated matrices and quadrature formulas
                % for triangles (1978)
                % Number of points: 16
                % Degree of exactness: 8
                
                a1=1/3;                     w1=0.144315607677787;
                a2=0.658861384496478;       w2=0.103217370534718;
                b2=0.170569307751761;       w3=0.032458497623198;
                a3=0.898905543365938;       w4=0.095091634267284;
                b3=0.050547228317031;       w5=0.027230314174435;
                a4=0.081414823414554;
                b4=0.459292588292723;
                a5=0.008394777409958;
                b5=0.263112829634638;
                c5=0.728492392955404;
                
                points=[a1 a1;
                        b2 b2;
                        b2 a2;
                        a2 b2;
                        b3 b3;
                        b3 a3;
                        a3 b3;
                        b4 b4;
                        b4 a4;
                        a4 b4;
                        a5 b5;
                        a5 c5;
                        b5 a5;
                        b5 c5;
                        c5 a5;
                        c5 b5];
                    
                weights=0.5*[w1 w2 w2 w2 w3 w3 w3 w4 w4 w4 w5 w5 w5 w5 w5 w5];
                
            case 9
                % Source: Laursen and Gellert, Some criteria for
                % numerically integrated matrices and quadrature formulas
                % for triangles (1978)
                % Number of points: 19
                % Degree of exactness: 9
                
                a1=1/3;                     w1=0.097135796282799;
                a2=0.020634961602524;       w2=0.031334700227139;
                b2=0.489682519198738;       w3=0.077827541004774;
                a3=0.125820817014126;       w4=0.079647738927210;
                b3=0.437089591492937;       w5=0.025577675658698;
                a4=0.623592928761934;       w6=0.043283539377289;
                b4=0.188203535619033;
                a5=0.910540973211094;
                b5=0.044729513394453;
                a6=0.036838412054736;
                b6=0.221962989160766;
                c6=0.741198598784498;
                
                points=[a1 a1;
                        b2 b2;
                        b2 a2;
                        a2 b2;
                        b3 b3;
                        b3 a3;
                        a3 b3;
                        b4 b4;
                        b4 a4;
                        a4 b4;
                        b5 b5;
                        b5 a5;
                        a5 b5;
                        a6 b6;
                        a6 c6;
                        b6 a6;
                        b6 c6;
                        c6 a6;
                        c6 b6];
              
             weights=0.5*[w1 w2 w2 w2 w3 w3 w3 w4 w4 w4 w5 w5 w5 w6 w6 w6 w6 w6 w6];       
                    
            case 10
                % Source: Laursen and Gellert, Some criteria for
                % numerically integrated matrices and quadrature formulas
                % for triangles (1978)
                % Number of points: 25
                % Degree of exactness: 10
                  
                a1=1/3;                     w1=0.081743329146286;
                a2=0.715677797886872;       w2=0.045957963604745;
                b2=0.142161101056564;       w3=0.013352968813150;
                a3=0.935889253566112;       w4=0.063904906396424;
                b3=0.032055373216944;       w5=0.034184648162959;
                a4=0.148132885783821;       w6=0.025297757707288;
                b4=0.321812995288835;
                c4=0.530054118927344;
                a5=0.029619889488730;
                b5=0.369146781827811;
                c5=0.601233328683459;
                a6=0.028367665339938;
                b6=0.163701733737182;
                c6=0.807930600922880;
                 
                points=[a1 a1;
                        b2 b2;
                        a2 b2;
                        b2 a2;
                        b3 b3;
                        a3 b3;
                        b3 a3;
                        a4 b4;
                        a4 c4;
                        b4 a4;
                        b4 c4;
                        c4 a4;
                        c4 b4;
                        a5 b5;
                        a5 c5;
                        b5 a5;
                        b5 c5;
                        c5 a5;
                        c5 b5;
                        a6 b6;
                        a6 c6;
                        b6 a6;
                        b6 c6;
                        c6 a6;
                        c6 b6];
                                    
                weights=0.5*[w1 w2 w2 w2 w3 w3 w3 w4 w4 w4 w4 w4 w4 w5 w5 w5 w5 w5 w5 w6 w6 w6 w6 w6 w6];           
                
            otherwise
                error('Unimplemented');
        end
        
    case 'boundary'
        
        [QuadRule] = getGaussLegendre(0, 1, id);
        
        % Quadrature nodes
        points=QuadRule.x;
        % Quadrature weights
        weights=QuadRule.w;
        
end


    function [QuadRule] = getGaussLegendre(a,b,n)
        
        % getGaussLegendre: 1D n-point Gauss-Legendre quadrature rule on the interval [a,b].
        % Quadrature rules obtained from Gauss-Legendre are of order 2*n-1.
        % INPUT:
        % a,b:          End-points of the interval
        % n:            Number of quadrature nodes in [a,b]
        % OUTPUT:
        % QuadRule:     Structure which contains the fields:
        %    x:         N-by-1 matrix specifying the abscissae of the quadrature rule.
        %    w:         N-by-1 matrix specifying the weights of the quadrature rule.
        
        if (n==1), x = 0; w = 2;
        else
            c = zeros(n-1,1);
            for i=1:(n-1), c(i)=i/sqrt(4*i*i-1); end
            J=diag(c,-1)+diag(c,1); [ev,ew]=eig(J);
            x=diag(ew); w=(2*(ev(1,:).*ev(1,:)))';
        end
        
        % The above rule computes a quadrature rule for the reference interval.
        % The reference interval is always [-1,1], and if [a,b] != [-1,1] then we
        % re-map the abscissas (x) and rescale the weights.
        
        % Midpoint
        xm = (b+a)/2;
        % Area of the requested interval / area of the reference interval
        xl = (b-a)/2;
        
        x = xm + xl*x;
        w = w * xl;
        
        QuadRule.w = w(:);
        QuadRule.x = x(:);
    end
end