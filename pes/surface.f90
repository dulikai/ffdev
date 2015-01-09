! pes for c3 grid by prof. Z. Lan
! 2014.12
      program main
      implicit none
      double precision ::   v_1, v_2, theta, rho

   

                       
      character (len=50) :: file1
      character (len=50) :: file2
      character (len=50) :: file3
      character (len=50) :: file4
      integer :: i, j, m, n, state
      integer :: n_i, n_j 
      double precision :: ome,K, F, G, q_x, q_y,  &
                          qx_min, qy_min, &
                          qx_max, qy_max, &
                          dqx, dqy
                

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Calculations of absorption spectra
!    By Z. Lan
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      ome = 1 
      K = 3.5
      F = 1
      G = 0.5 


      qx_min = -5
      qy_min = -5
      qx_max =  5
      qy_max =  5
      n_i = 200
      n_j = 200
      dqx = ( qx_max - qx_min) / n_i
      dqy = ( qy_max - qy_min) / n_j
        

      v_1 = 0
      v_2 = 0
      theta = 0 
      rho = 0
      do i = 1 , n_i
      do j = 1 , n_j
         q_x = qx_min + (i - 1) * dqx
         q_y = qy_min + (j - 1) * dqy
         rho = (q_x**2 + q_y**2) **0.5
         
         if ( (q_x  .ge. 0 ) .AND.     &
              (q_y  .ge. 0 )           & 
            )  then
             theta = atan(q_y/q_x) 
         endif
         if ( (q_x  .lt. 0 ) .AND.     &
              (q_y  .ge. 0 )           &
            ) then
             theta = 3.1415 + atan(q_y/q_x)
         endif
          if ( (q_x  .lt. 0 ) .AND.     &
              (q_y  .lt. 0 )           &
            ) then
             theta = 3.1415 + atan(q_y/q_x)
         endif
          if ( (q_x  .ge. 0 ) .AND.     &
              (q_y   .lt. 0 )           &
            ) then
             theta = 6.283 + atan(q_y/q_x)
         endif

         v_1      =  0.5 * ome * rho **2    &
                     - 0.5* ( 4*K**2*rho**2      &
                              + 2*K*G*rho**3*cos(3*theta)  &
                              + G**2*rho**4 )**0.5        
         v_2      = 0.5 * ome * rho **2    &
                     + 0.5* ( 4*K**2*rho**2      &
                               + 2*K*G*rho**3*cos(3*theta) &
                              + G**2*rho**4 )**0.5

     
         write (*,9999) q_x, q_y, rho, theta/3.14*180, v_1, v_2
        
9999     format (10(f20.10,1x))

      enddo 
!        write (*,*)
      enddo



       end
