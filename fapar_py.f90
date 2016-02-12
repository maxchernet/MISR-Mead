program fapar

end program fapar


! Uncertainty

subroutine fapar_sd(f_name_in, f_name_out)
        PARAMETER (PI=3.141592653589793)
        INTEGER :: LAD, NV, stat, sample_n, iend, year
        REAL :: BRF(1000)
        REAL :: THETA_V,PHI_V,THETA_I,PHI_I
        REAL :: XLAI,XHC,RPL,XRL,XTL,XRS
        REAL :: test_nv
        real :: teta_s,phi_s
        character*6 :: day_str, year_str
        character*100 :: f_name_in, f_name_out

!        write(day_str,'(I6)')i
        iend = len(day_str)
!        f_name_in = '/home/max/s3vt_ng/output_misr_etm/misr_etm/output_fapar_sd/fapar_in_'&
!        //year_str(index(year_str,' ', back=.true.)+1:iend)//'_'&
!        &//day_str(index(day_str,' ', back=.true.)+1:iend)//'_ang7_Ne1.dat'
        OPEN(10, file=f_name_in, STATUS='OLD')
!        f_name_out = '/home/max/s3vt_ng/output_misr_etm/misr_etm/output_fapar_nadim_sd/fapar_out_Ne1_'&
!        //year_str(index(year_str,' ', back=.true.)+1:iend)//'_'&
!        &//day_str(index(day_str,' ', back=.true.)+1:iend)//'_ang7.dat'
        open(11, file=f_name_out)
        sample_n=1
        do
                READ(10, *, iostat=stat) theta_i, ph_i, nv, theta_v, phi_v, lad, xrs, xhc, xlai, rpl, xrl, xtl
                if (stat < 0) exit
                call NADIMBRF(THETA_I, PHI_I, NV, THETA_V, PHI_V, LAD, XRS, XHC, XLAI, RPL, XRL, XTL, BRF)
                call energie(theta_i, phi_i, fpar, albedo_sys, trans_totale)
                ! print*, theta_i, ph_i, nv, theta_v, phi_v, lad, xrs, xhc, xlai, rpl, xrl, xtl
                ! print*, '\n'
                print*,'fpar=', sample_n, fpar
                sample_n = sample_n + 1
                WRITE(11,'(F8.3)') fpar
        end do
        close(10)
        close(11)

end subroutine fapar_sd


! Mean

subroutine fapar_mean(f_name_in, f_name_out)
        PARAMETER (PI=3.141592653589793)
        INTEGER :: LAD, NV, stat, sample_n, iend, year
        REAL :: BRF(1000)
        REAL :: THETA_V,PHI_V,THETA_I,PHI_I
        REAL :: XLAI,XHC,RPL,XRL,XTL,XRS
        REAL :: test_nv
        real :: teta_s,phi_s
        character*6 :: day_str, year_str
        character*100 :: f_name_in, f_name_out

!        write(year_str,'(I6)')year
!        iend = len(year_str)
!        f_name_in = '/home/max/s3vt_ng/output_misr_etm/misr_etm/output_fapar_mean/fapar_in_Ne1_'&
!        //year_str(index(year_str,' ', back=.true.)+1:iend)//'_mean.dat'
        OPEN(10, file=f_name_in, STATUS='OLD')
!        f_name_out = '/home/max/s3vt_ng/output_misr_etm/misr_etm/output_fapar_nadim_mean/fapar_out_Ne1_'&
!        //year_str(index(year_str,' ', back=.true.)+1:iend)//'_mean.dat'
        open(11, file=f_name_out)
        day=1
        !print*, f_name_in, '\n'
        do
                READ(10, *, iostat=stat) theta_i, ph_i, nv, theta_v, phi_v, lad, xrs, xhc, xlai, rpl, xrl, xtl
! Check end of file
                if (stat < 0) exit
                ! print*, theta_i, ph_i, nv, theta_v, phi_v, lad, xrs, xhc, xlai, rpl, xrl, xtl
                !print*, '\n'
                call NADIMBRF(THETA_I,PHI_I,NV,THETA_V,PHI_V,LAD,XRS,XHC,XLAI,RPL,XRL,XTL,BRF)
                call energie(theta_i,phi_i,fpar,albedo_sys,trans_totale)
                print*,'mean fpar=', day, fpar
                day = day + 1
                WRITE(11,'(F8.3)') fpar
        end do
        close(10)
        close(11)
        return
end subroutine fapar_mean
