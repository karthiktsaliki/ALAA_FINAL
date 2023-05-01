function B = Bidiag_Francis_Step( B )
    %Francis_Step Perform one Francis Implicit QR Step
    B;
    [ m, n ] = size( B );
    % if m <= 2 
    %     return
    % end
    
    % Introduce the bulge
    % Compute the first Givens' rotation


    Gt0 = Givens_rotation( [ B(1,1)^2 - (B(m-1,m)^2+B(m,m)^2) 
                             B( 1,2 )*B(1,1) ]);
    Gt0;
   
    B(1:2,1:2) = B(1:2,1:2)*Gt0;

    for i=1:m-1
      F0 =  Givens_rotation([ B(i,i)
                               B(i+1,i)]);
      if i+2 <= m
          B(i:i+1,i:i+2) = F0'*B(i:i+1,i:i+2);
      else
          B(i:i+1,i:i+1) = F0'*B(i:i+1,i:i+1);

      end
      if (i<m-1)
          Gt1 = Givens_rotation([B(i,i+1)
                                 B(i,i+2)]);
          B(i:i+2,i+1:i+2)=B(i:i+2,i+1:i+2)*Gt1;
      end
     
    end
end
