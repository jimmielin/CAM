
      module mo_nln_matrix

      use shr_kind_mod, only : r8 => shr_kind_r8

      use chem_mods, only: veclen
      private
      public :: nlnmat

      contains

      subroutine     nlnmat01( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,638) = -(rxt(k,356)*y(k,223))
         mat(k,1637) = -rxt(k,356)*y(k,1)

         mat(k,1806) = rxt(k,359)*y(k,195)
         mat(k,917) = rxt(k,359)*y(k,130)

         mat(k,672) = -(rxt(k,360)*y(k,223))
         mat(k,1640) = -rxt(k,360)*y(k,2)

         mat(k,918) = rxt(k,357)*y(k,209)
         mat(k,1959) = rxt(k,357)*y(k,195)




         mat(k,976) = -(rxt(k,439)*y(k,132) + rxt(k,440)*y(k,140) + rxt(k,441) &
                      *y(k,223))
         mat(k,1728) = -rxt(k,439)*y(k,6)
         mat(k,2132) = -rxt(k,440)*y(k,6)
         mat(k,1668) = -rxt(k,441)*y(k,6)

         mat(k,168) = -(rxt(k,398)*y(k,223))
         mat(k,1568) = -rxt(k,398)*y(k,7)

         mat(k,397) = -(rxt(k,401)*y(k,223))
         mat(k,1605) = -rxt(k,401)*y(k,8)

         mat(k,495) = rxt(k,399)*y(k,209)
         mat(k,1938) = rxt(k,399)*y(k,197)


         mat(k,169) = .120_r8*rxt(k,398)*y(k,223)
         mat(k,1569) = .120_r8*rxt(k,398)*y(k,7)


         mat(k,970) = .100_r8*rxt(k,440)*y(k,140)
         mat(k,895) = .100_r8*rxt(k,443)*y(k,140)
         mat(k,2117) = .100_r8*rxt(k,440)*y(k,6) + .100_r8*rxt(k,443)*y(k,116)


         mat(k,1793) = .500_r8*rxt(k,400)*y(k,197) + .200_r8*rxt(k,427)*y(k,229)  &
                      + .060_r8*rxt(k,433)*y(k,232)
         mat(k,496) = .500_r8*rxt(k,400)*y(k,130)
         mat(k,734) = .200_r8*rxt(k,427)*y(k,130)
         mat(k,750) = .060_r8*rxt(k,433)*y(k,130)


         mat(k,1787) = .200_r8*rxt(k,427)*y(k,229) + .200_r8*rxt(k,433)*y(k,232)
         mat(k,733) = .200_r8*rxt(k,427)*y(k,130)
         mat(k,748) = .200_r8*rxt(k,433)*y(k,130)


         mat(k,1803) = .200_r8*rxt(k,427)*y(k,229) + .150_r8*rxt(k,433)*y(k,232)
         mat(k,735) = .200_r8*rxt(k,427)*y(k,130)
         mat(k,751) = .150_r8*rxt(k,433)*y(k,130)


         mat(k,1789) = .210_r8*rxt(k,433)*y(k,232)
         mat(k,749) = .210_r8*rxt(k,433)*y(k,130)

         mat(k,250) = -(rxt(k,361)*y(k,223))
         mat(k,1583) = -rxt(k,361)*y(k,15)

         mat(k,969) = .050_r8*rxt(k,440)*y(k,140)
         mat(k,894) = .050_r8*rxt(k,443)*y(k,140)
         mat(k,2116) = .050_r8*rxt(k,440)*y(k,6) + .050_r8*rxt(k,443)*y(k,116)

         mat(k,355) = -(rxt(k,327)*y(k,132) + rxt(k,328)*y(k,223))
         mat(k,1718) = -rxt(k,327)*y(k,16)
         mat(k,1599) = -rxt(k,328)*y(k,16)

         mat(k,1423) = -(rxt(k,211)*y(k,42) + rxt(k,212)*y(k,209) + rxt(k,213) &
                      *y(k,140))
         mat(k,2027) = -rxt(k,211)*y(k,17)
         mat(k,2003) = -rxt(k,212)*y(k,17)
         mat(k,2154) = -rxt(k,213)*y(k,17)

         mat(k,1490) = 4.000_r8*rxt(k,214)*y(k,19) + (rxt(k,215)+rxt(k,216))*y(k,59)  &
                      + rxt(k,219)*y(k,130) + rxt(k,222)*y(k,139) + rxt(k,469) &
                      *y(k,156) + rxt(k,223)*y(k,223)
         mat(k,149) = rxt(k,201)*y(k,222)
         mat(k,155) = rxt(k,227)*y(k,222)
         mat(k,474) = 2.000_r8*rxt(k,238)*y(k,56) + 2.000_r8*rxt(k,250)*y(k,222)  &
                      + 2.000_r8*rxt(k,239)*y(k,223)
         mat(k,592) = rxt(k,240)*y(k,56) + rxt(k,251)*y(k,222) + rxt(k,241)*y(k,223)
         mat(k,454) = 3.000_r8*rxt(k,245)*y(k,56) + 3.000_r8*rxt(k,228)*y(k,222)  &
                      + 3.000_r8*rxt(k,246)*y(k,223)
         mat(k,2092) = 2.000_r8*rxt(k,238)*y(k,41) + rxt(k,240)*y(k,43)  &
                      + 3.000_r8*rxt(k,245)*y(k,55)
         mat(k,2054) = (rxt(k,215)+rxt(k,216))*y(k,19)
         mat(k,113) = 2.000_r8*rxt(k,229)*y(k,222)
         mat(k,816) = rxt(k,224)*y(k,139) + rxt(k,230)*y(k,222) + rxt(k,225)*y(k,223)
         mat(k,1845) = rxt(k,219)*y(k,19)
         mat(k,2252) = rxt(k,222)*y(k,19) + rxt(k,224)*y(k,87)
         mat(k,1238) = rxt(k,469)*y(k,19)
         mat(k,1530) = rxt(k,201)*y(k,34) + rxt(k,227)*y(k,35) + 2.000_r8*rxt(k,250) &
                      *y(k,41) + rxt(k,251)*y(k,43) + 3.000_r8*rxt(k,228)*y(k,55)  &
                      + 2.000_r8*rxt(k,229)*y(k,84) + rxt(k,230)*y(k,87)
         mat(k,1694) = rxt(k,223)*y(k,19) + 2.000_r8*rxt(k,239)*y(k,41) + rxt(k,241) &
                      *y(k,43) + 3.000_r8*rxt(k,246)*y(k,55) + rxt(k,225)*y(k,87)


         mat(k,1484) = rxt(k,217)*y(k,59)
         mat(k,2048) = rxt(k,217)*y(k,19)
         mat(k,2174) = (rxt(k,537)+rxt(k,542))*y(k,97)
         mat(k,783) = (rxt(k,537)+rxt(k,542))*y(k,91)

         mat(k,1493) = -(4._r8*rxt(k,214)*y(k,19) + (rxt(k,215) + rxt(k,216) + rxt(k,217) &
                      ) * y(k,59) + rxt(k,218)*y(k,209) + rxt(k,219)*y(k,130) &
                      + rxt(k,220)*y(k,131) + rxt(k,222)*y(k,139) + rxt(k,223) &
                      *y(k,223) + rxt(k,469)*y(k,156))
         mat(k,2057) = -(rxt(k,215) + rxt(k,216) + rxt(k,217)) * y(k,19)
         mat(k,2007) = -rxt(k,218)*y(k,19)
         mat(k,1849) = -rxt(k,219)*y(k,19)
         mat(k,2225) = -rxt(k,220)*y(k,19)
         mat(k,2256) = -rxt(k,222)*y(k,19)
         mat(k,1698) = -rxt(k,223)*y(k,19)
         mat(k,1240) = -rxt(k,469)*y(k,19)

         mat(k,1425) = rxt(k,213)*y(k,140)
         mat(k,566) = rxt(k,221)*y(k,139)
         mat(k,818) = rxt(k,231)*y(k,222)
         mat(k,786) = rxt(k,226)*y(k,139)
         mat(k,2256) = mat(k,2256) + rxt(k,221)*y(k,20) + rxt(k,226)*y(k,97)
         mat(k,2158) = rxt(k,213)*y(k,17)
         mat(k,1534) = rxt(k,231)*y(k,87)

         mat(k,562) = -(rxt(k,221)*y(k,139))
         mat(k,2242) = -rxt(k,221)*y(k,20)

         mat(k,1486) = rxt(k,220)*y(k,131)
         mat(k,2204) = rxt(k,220)*y(k,19)


         mat(k,244) = -(rxt(k,402)*y(k,223))
         mat(k,1581) = -rxt(k,402)*y(k,22)

         mat(k,1784) = rxt(k,405)*y(k,199)
         mat(k,433) = rxt(k,405)*y(k,130)

         mat(k,327) = -(rxt(k,404)*y(k,223))
         mat(k,1595) = -rxt(k,404)*y(k,23)

         mat(k,434) = rxt(k,403)*y(k,209)
         mat(k,1933) = rxt(k,403)*y(k,199)

         mat(k,285) = -(rxt(k,276)*y(k,56) + rxt(k,277)*y(k,223))
         mat(k,2073) = -rxt(k,276)*y(k,24)
         mat(k,1589) = -rxt(k,277)*y(k,24)

         mat(k,554) = -(rxt(k,278)*y(k,56) + rxt(k,279)*y(k,140) + rxt(k,304)*y(k,223))
         mat(k,2078) = -rxt(k,278)*y(k,25)
         mat(k,2120) = -rxt(k,279)*y(k,25)
         mat(k,1627) = -rxt(k,304)*y(k,25)

         mat(k,262) = -(rxt(k,284)*y(k,223))
         mat(k,1586) = -rxt(k,284)*y(k,26)

         mat(k,823) = .800_r8*rxt(k,280)*y(k,200) + .200_r8*rxt(k,281)*y(k,204)
         mat(k,1864) = .200_r8*rxt(k,281)*y(k,200)

         mat(k,345) = -(rxt(k,285)*y(k,223))
         mat(k,1598) = -rxt(k,285)*y(k,27)

         mat(k,824) = rxt(k,282)*y(k,209)
         mat(k,1935) = rxt(k,282)*y(k,200)

         mat(k,294) = -(rxt(k,286)*y(k,56) + rxt(k,287)*y(k,223))
         mat(k,2074) = -rxt(k,286)*y(k,28)
         mat(k,1590) = -rxt(k,287)*y(k,28)

         mat(k,1009) = -(rxt(k,307)*y(k,132) + rxt(k,308)*y(k,140) + rxt(k,325) &
                      *y(k,223))
         mat(k,1730) = -rxt(k,307)*y(k,29)
         mat(k,2134) = -rxt(k,308)*y(k,29)
         mat(k,1670) = -rxt(k,325)*y(k,29)

         mat(k,855) = .130_r8*rxt(k,385)*y(k,140)
         mat(k,2134) = mat(k,2134) + .130_r8*rxt(k,385)*y(k,104)

         mat(k,409) = -(rxt(k,312)*y(k,223))
         mat(k,1607) = -rxt(k,312)*y(k,30)

         mat(k,792) = rxt(k,310)*y(k,209)
         mat(k,1940) = rxt(k,310)*y(k,201)

         mat(k,115) = -(rxt(k,313)*y(k,223))
         mat(k,1565) = -rxt(k,313)*y(k,31)

         mat(k,271) = -(rxt(k,408)*y(k,223))
         mat(k,1588) = -rxt(k,408)*y(k,32)

         mat(k,629) = rxt(k,406)*y(k,209)
         mat(k,1929) = rxt(k,406)*y(k,202)

         mat(k,105) = -(rxt(k,200)*y(k,222))
         mat(k,1508) = -rxt(k,200)*y(k,33)

         mat(k,147) = -(rxt(k,201)*y(k,222))
         mat(k,1513) = -rxt(k,201)*y(k,34)

         mat(k,152) = -(rxt(k,227)*y(k,222))
         mat(k,1514) = -rxt(k,227)*y(k,35)

         mat(k,119) = -(rxt(k,202)*y(k,222))
         mat(k,1510) = -rxt(k,202)*y(k,36)

         mat(k,157) = -(rxt(k,203)*y(k,222))
         mat(k,1515) = -rxt(k,203)*y(k,37)

         mat(k,123) = -(rxt(k,204)*y(k,222))
         mat(k,1511) = -rxt(k,204)*y(k,38)

         mat(k,162) = -(rxt(k,205)*y(k,222))
         mat(k,1516) = -rxt(k,205)*y(k,39)

         mat(k,127) = -(rxt(k,206)*y(k,222))
         mat(k,1512) = -rxt(k,206)*y(k,40)

         mat(k,473) = -(rxt(k,238)*y(k,56) + rxt(k,239)*y(k,223) + rxt(k,250)*y(k,222))
         mat(k,2077) = -rxt(k,238)*y(k,41)
         mat(k,1616) = -rxt(k,239)*y(k,41)
         mat(k,1525) = -rxt(k,250)*y(k,41)

         mat(k,2038) = -(rxt(k,175)*y(k,56) + rxt(k,211)*y(k,17) + rxt(k,255)*y(k,209) &
                      + rxt(k,256)*y(k,132) + rxt(k,257)*y(k,139) + rxt(k,258) &
                      *y(k,223))
         mat(k,2103) = -rxt(k,175)*y(k,42)
         mat(k,1429) = -rxt(k,211)*y(k,42)
         mat(k,2014) = -rxt(k,255)*y(k,42)
         mat(k,1764) = -rxt(k,256)*y(k,42)
         mat(k,2263) = -rxt(k,257)*y(k,42)
         mat(k,1705) = -rxt(k,258)*y(k,42)

         mat(k,646) = .400_r8*rxt(k,356)*y(k,223)
         mat(k,991) = .340_r8*rxt(k,440)*y(k,140)
         mat(k,361) = .500_r8*rxt(k,327)*y(k,132)
         mat(k,560) = rxt(k,279)*y(k,140)
         mat(k,1021) = .500_r8*rxt(k,308)*y(k,140)
         mat(k,621) = .500_r8*rxt(k,296)*y(k,223)
         mat(k,806) = rxt(k,263)*y(k,223)
         mat(k,451) = .300_r8*rxt(k,264)*y(k,223)
         mat(k,1446) = (rxt(k,272)+rxt(k,273))*y(k,222)
         mat(k,2064) = rxt(k,182)*y(k,204)
         mat(k,1057) = .800_r8*rxt(k,301)*y(k,223)
         mat(k,867) = .910_r8*rxt(k,385)*y(k,140)
         mat(k,585) = .300_r8*rxt(k,376)*y(k,223)
         mat(k,1209) = .800_r8*rxt(k,380)*y(k,204)
         mat(k,1221) = .120_r8*rxt(k,338)*y(k,140)
         mat(k,614) = .500_r8*rxt(k,351)*y(k,223)
         mat(k,913) = .340_r8*rxt(k,443)*y(k,140)
         mat(k,1350) = .600_r8*rxt(k,352)*y(k,140)
         mat(k,1856) = .100_r8*rxt(k,358)*y(k,195) + rxt(k,262)*y(k,204)  &
                      + .500_r8*rxt(k,329)*y(k,206) + .500_r8*rxt(k,298)*y(k,208)  &
                      + .920_r8*rxt(k,368)*y(k,211) + .250_r8*rxt(k,336)*y(k,215)  &
                      + rxt(k,345)*y(k,217) + rxt(k,319)*y(k,225) + rxt(k,323) &
                      *y(k,226) + .340_r8*rxt(k,452)*y(k,227) + .320_r8*rxt(k,457) &
                      *y(k,228) + .250_r8*rxt(k,393)*y(k,231)
         mat(k,1764) = mat(k,1764) + .500_r8*rxt(k,327)*y(k,16) + rxt(k,369)*y(k,211)  &
                      + .250_r8*rxt(k,335)*y(k,215) + rxt(k,346)*y(k,217)
         mat(k,2165) = .340_r8*rxt(k,440)*y(k,6) + rxt(k,279)*y(k,25)  &
                      + .500_r8*rxt(k,308)*y(k,29) + .910_r8*rxt(k,385)*y(k,104)  &
                      + .120_r8*rxt(k,338)*y(k,111) + .340_r8*rxt(k,443)*y(k,116)  &
                      + .600_r8*rxt(k,352)*y(k,117)
         mat(k,536) = rxt(k,303)*y(k,223)
         mat(k,1084) = .680_r8*rxt(k,461)*y(k,223)
         mat(k,929) = .100_r8*rxt(k,358)*y(k,130)
         mat(k,832) = .700_r8*rxt(k,281)*y(k,204)
         mat(k,800) = rxt(k,309)*y(k,204)
         mat(k,1402) = rxt(k,292)*y(k,204) + rxt(k,365)*y(k,211) + .250_r8*rxt(k,332) &
                      *y(k,215) + rxt(k,341)*y(k,217) + .250_r8*rxt(k,390)*y(k,231)
         mat(k,1907) = rxt(k,182)*y(k,59) + .800_r8*rxt(k,380)*y(k,107) + rxt(k,262) &
                      *y(k,130) + .700_r8*rxt(k,281)*y(k,200) + rxt(k,309)*y(k,201)  &
                      + rxt(k,292)*y(k,203) + (4.000_r8*rxt(k,259)+2.000_r8*rxt(k,260)) &
                      *y(k,204) + 1.500_r8*rxt(k,366)*y(k,211) + .750_r8*rxt(k,371) &
                      *y(k,212) + .880_r8*rxt(k,333)*y(k,215) + 2.000_r8*rxt(k,342) &
                      *y(k,217) + .750_r8*rxt(k,445)*y(k,221) + .800_r8*rxt(k,321) &
                      *y(k,226) + .930_r8*rxt(k,450)*y(k,227) + .950_r8*rxt(k,455) &
                      *y(k,228) + .800_r8*rxt(k,391)*y(k,231)
         mat(k,576) = .500_r8*rxt(k,329)*y(k,130)
         mat(k,725) = .500_r8*rxt(k,298)*y(k,130)
         mat(k,2014) = mat(k,2014) + .450_r8*rxt(k,343)*y(k,217) + .150_r8*rxt(k,322) &
                      *y(k,226)
         mat(k,1273) = .920_r8*rxt(k,368)*y(k,130) + rxt(k,369)*y(k,132) + rxt(k,365) &
                      *y(k,203) + 1.500_r8*rxt(k,366)*y(k,204)
         mat(k,1306) = .750_r8*rxt(k,371)*y(k,204)
         mat(k,1328) = .250_r8*rxt(k,336)*y(k,130) + .250_r8*rxt(k,335)*y(k,132)  &
                      + .250_r8*rxt(k,332)*y(k,203) + .880_r8*rxt(k,333)*y(k,204)
         mat(k,1370) = rxt(k,345)*y(k,130) + rxt(k,346)*y(k,132) + rxt(k,341)*y(k,203)  &
                      + 2.000_r8*rxt(k,342)*y(k,204) + .450_r8*rxt(k,343)*y(k,209)  &
                      + 4.000_r8*rxt(k,344)*y(k,217)
         mat(k,1073) = .750_r8*rxt(k,445)*y(k,204)
         mat(k,1541) = (rxt(k,272)+rxt(k,273))*y(k,54)
         mat(k,1705) = mat(k,1705) + .400_r8*rxt(k,356)*y(k,1) + .500_r8*rxt(k,296) &
                      *y(k,51) + rxt(k,263)*y(k,52) + .300_r8*rxt(k,264)*y(k,53)  &
                      + .800_r8*rxt(k,301)*y(k,80) + .300_r8*rxt(k,376)*y(k,105)  &
                      + .500_r8*rxt(k,351)*y(k,115) + rxt(k,303)*y(k,145)  &
                      + .680_r8*rxt(k,461)*y(k,184)
         mat(k,780) = rxt(k,319)*y(k,130)
         mat(k,1168) = rxt(k,323)*y(k,130) + .800_r8*rxt(k,321)*y(k,204)  &
                      + .150_r8*rxt(k,322)*y(k,209)
         mat(k,1132) = .340_r8*rxt(k,452)*y(k,130) + .930_r8*rxt(k,450)*y(k,204)
         mat(k,1154) = .320_r8*rxt(k,457)*y(k,130) + .950_r8*rxt(k,455)*y(k,204)
         mat(k,1186) = .250_r8*rxt(k,393)*y(k,130) + .250_r8*rxt(k,390)*y(k,203)  &
                      + .800_r8*rxt(k,391)*y(k,204)

      end do

      end subroutine     nlnmat01

      subroutine     nlnmat02( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,591) = -(rxt(k,240)*y(k,56) + rxt(k,241)*y(k,223) + rxt(k,251)*y(k,222))
         mat(k,2079) = -rxt(k,240)*y(k,43)
         mat(k,1631) = -rxt(k,241)*y(k,43)
         mat(k,1526) = -rxt(k,251)*y(k,43)

         mat(k,131) = -(rxt(k,242)*y(k,223))
         mat(k,1566) = -rxt(k,242)*y(k,44)

         mat(k,1042) = -(rxt(k,288)*y(k,132) + rxt(k,289)*y(k,223))
         mat(k,1732) = -rxt(k,288)*y(k,45)
         mat(k,1672) = -rxt(k,289)*y(k,45)

         mat(k,642) = .800_r8*rxt(k,356)*y(k,223)
         mat(k,358) = rxt(k,327)*y(k,132)
         mat(k,263) = rxt(k,284)*y(k,223)
         mat(k,347) = .500_r8*rxt(k,285)*y(k,223)
         mat(k,1010) = .500_r8*rxt(k,308)*y(k,140)
         mat(k,1335) = .100_r8*rxt(k,352)*y(k,140)
         mat(k,1825) = .400_r8*rxt(k,358)*y(k,195) + rxt(k,283)*y(k,200)  &
                      + .270_r8*rxt(k,311)*y(k,201) + rxt(k,329)*y(k,206) + rxt(k,348) &
                      *y(k,219) + rxt(k,319)*y(k,225)
         mat(k,1732) = mat(k,1732) + rxt(k,327)*y(k,16)
         mat(k,2135) = .500_r8*rxt(k,308)*y(k,29) + .100_r8*rxt(k,352)*y(k,117)
         mat(k,923) = .400_r8*rxt(k,358)*y(k,130)
         mat(k,827) = rxt(k,283)*y(k,130) + 3.200_r8*rxt(k,280)*y(k,200)  &
                      + .800_r8*rxt(k,281)*y(k,204)
         mat(k,795) = .270_r8*rxt(k,311)*y(k,130)
         mat(k,1879) = .800_r8*rxt(k,281)*y(k,200)
         mat(k,572) = rxt(k,329)*y(k,130)
         mat(k,1983) = .200_r8*rxt(k,347)*y(k,219)
         mat(k,684) = rxt(k,348)*y(k,130) + .200_r8*rxt(k,347)*y(k,209)
         mat(k,1672) = mat(k,1672) + .800_r8*rxt(k,356)*y(k,1) + rxt(k,284)*y(k,26)  &
                      + .500_r8*rxt(k,285)*y(k,27)
         mat(k,775) = rxt(k,319)*y(k,130)

         mat(k,371) = -(rxt(k,243)*y(k,56) + rxt(k,244)*y(k,223))
         mat(k,2075) = -rxt(k,243)*y(k,46)
         mat(k,1601) = -rxt(k,244)*y(k,46)

         mat(k,108) = -(rxt(k,290)*y(k,223))
         mat(k,1564) = -rxt(k,290)*y(k,47)

         mat(k,950) = -(rxt(k,326)*y(k,223))
         mat(k,1666) = -rxt(k,326)*y(k,48)

         mat(k,641) = .800_r8*rxt(k,356)*y(k,223)
         mat(k,974) = .520_r8*rxt(k,440)*y(k,140)
         mat(k,357) = .500_r8*rxt(k,327)*y(k,132)
         mat(k,900) = .520_r8*rxt(k,443)*y(k,140)
         mat(k,1821) = .250_r8*rxt(k,358)*y(k,195) + .820_r8*rxt(k,311)*y(k,201)  &
                      + .500_r8*rxt(k,329)*y(k,206) + .270_r8*rxt(k,452)*y(k,227)  &
                      + .040_r8*rxt(k,457)*y(k,228)
         mat(k,1726) = .500_r8*rxt(k,327)*y(k,16)
         mat(k,2130) = .520_r8*rxt(k,440)*y(k,6) + .520_r8*rxt(k,443)*y(k,116)
         mat(k,1076) = .500_r8*rxt(k,461)*y(k,223)
         mat(k,922) = .250_r8*rxt(k,358)*y(k,130)
         mat(k,794) = .820_r8*rxt(k,311)*y(k,130) + .820_r8*rxt(k,309)*y(k,204)
         mat(k,1875) = .820_r8*rxt(k,309)*y(k,201) + .150_r8*rxt(k,450)*y(k,227)  &
                      + .025_r8*rxt(k,455)*y(k,228)
         mat(k,571) = .500_r8*rxt(k,329)*y(k,130)
         mat(k,1666) = mat(k,1666) + .800_r8*rxt(k,356)*y(k,1) + .500_r8*rxt(k,461) &
                      *y(k,184)
         mat(k,1120) = .270_r8*rxt(k,452)*y(k,130) + .150_r8*rxt(k,450)*y(k,204)
         mat(k,1139) = .040_r8*rxt(k,457)*y(k,130) + .025_r8*rxt(k,455)*y(k,204)

         mat(k,1226) = -(rxt(k,314)*y(k,132) + rxt(k,315)*y(k,223))
         mat(k,1745) = -rxt(k,314)*y(k,49)
         mat(k,1685) = -rxt(k,315)*y(k,49)

         mat(k,1111) = rxt(k,316)*y(k,223)
         mat(k,1215) = .880_r8*rxt(k,338)*y(k,140)
         mat(k,1338) = .500_r8*rxt(k,352)*y(k,140)
         mat(k,1838) = .170_r8*rxt(k,411)*y(k,205) + .050_r8*rxt(k,374)*y(k,212)  &
                      + .250_r8*rxt(k,336)*y(k,215) + .170_r8*rxt(k,417)*y(k,218)  &
                      + .400_r8*rxt(k,427)*y(k,229) + .250_r8*rxt(k,393)*y(k,231)  &
                      + .540_r8*rxt(k,433)*y(k,232) + .510_r8*rxt(k,436)*y(k,234)
         mat(k,1745) = mat(k,1745) + .050_r8*rxt(k,375)*y(k,212) + .250_r8*rxt(k,335) &
                      *y(k,215) + .250_r8*rxt(k,394)*y(k,231)
         mat(k,838) = rxt(k,317)*y(k,223)
         mat(k,2146) = .880_r8*rxt(k,338)*y(k,111) + .500_r8*rxt(k,352)*y(k,117)
         mat(k,1388) = .250_r8*rxt(k,332)*y(k,215) + .250_r8*rxt(k,390)*y(k,231)
         mat(k,1891) = .240_r8*rxt(k,333)*y(k,215) + .500_r8*rxt(k,321)*y(k,226)  &
                      + .100_r8*rxt(k,391)*y(k,231)
         mat(k,767) = .170_r8*rxt(k,411)*y(k,130) + .070_r8*rxt(k,410)*y(k,209)
         mat(k,1995) = .070_r8*rxt(k,410)*y(k,205) + .070_r8*rxt(k,416)*y(k,218)
         mat(k,1294) = .050_r8*rxt(k,374)*y(k,130) + .050_r8*rxt(k,375)*y(k,132)
         mat(k,1318) = .250_r8*rxt(k,336)*y(k,130) + .250_r8*rxt(k,335)*y(k,132)  &
                      + .250_r8*rxt(k,332)*y(k,203) + .240_r8*rxt(k,333)*y(k,204)
         mat(k,872) = .170_r8*rxt(k,417)*y(k,130) + .070_r8*rxt(k,416)*y(k,209)
         mat(k,1685) = mat(k,1685) + rxt(k,316)*y(k,101) + rxt(k,317)*y(k,133)
         mat(k,1162) = .500_r8*rxt(k,321)*y(k,204)
         mat(k,743) = .400_r8*rxt(k,427)*y(k,130)
         mat(k,1179) = .250_r8*rxt(k,393)*y(k,130) + .250_r8*rxt(k,394)*y(k,132)  &
                      + .250_r8*rxt(k,390)*y(k,203) + .100_r8*rxt(k,391)*y(k,204)
         mat(k,759) = .540_r8*rxt(k,433)*y(k,130)
         mat(k,507) = .510_r8*rxt(k,436)*y(k,130)

         mat(k,690) = -(rxt(k,295)*y(k,223))
         mat(k,1642) = -rxt(k,295)*y(k,50)

         mat(k,1004) = .120_r8*rxt(k,308)*y(k,140)
         mat(k,2122) = .120_r8*rxt(k,308)*y(k,29)
         mat(k,1378) = .100_r8*rxt(k,292)*y(k,204) + .150_r8*rxt(k,293)*y(k,209)
         mat(k,1869) = .100_r8*rxt(k,292)*y(k,203)
         mat(k,1961) = .150_r8*rxt(k,293)*y(k,203) + .150_r8*rxt(k,343)*y(k,217)
         mat(k,1357) = .150_r8*rxt(k,343)*y(k,209)

         mat(k,616) = -(rxt(k,296)*y(k,223))
         mat(k,1634) = -rxt(k,296)*y(k,51)

         mat(k,1377) = .400_r8*rxt(k,293)*y(k,209)
         mat(k,1955) = .400_r8*rxt(k,293)*y(k,203) + .400_r8*rxt(k,343)*y(k,217)
         mat(k,1356) = .400_r8*rxt(k,343)*y(k,209)

         mat(k,803) = -(rxt(k,263)*y(k,223))
         mat(k,1652) = -rxt(k,263)*y(k,52)

         mat(k,1192) = .200_r8*rxt(k,380)*y(k,204)
         mat(k,825) = .300_r8*rxt(k,281)*y(k,204)
         mat(k,1871) = .200_r8*rxt(k,380)*y(k,107) + .300_r8*rxt(k,281)*y(k,200)  &
                      + 2.000_r8*rxt(k,260)*y(k,204) + .250_r8*rxt(k,366)*y(k,211)  &
                      + .250_r8*rxt(k,371)*y(k,212) + .250_r8*rxt(k,333)*y(k,215)  &
                      + .250_r8*rxt(k,445)*y(k,221) + .500_r8*rxt(k,321)*y(k,226)  &
                      + .250_r8*rxt(k,450)*y(k,227) + .250_r8*rxt(k,455)*y(k,228)  &
                      + .300_r8*rxt(k,391)*y(k,231)
         mat(k,1252) = .250_r8*rxt(k,366)*y(k,204)
         mat(k,1283) = .250_r8*rxt(k,371)*y(k,204)
         mat(k,1312) = .250_r8*rxt(k,333)*y(k,204)
         mat(k,1061) = .250_r8*rxt(k,445)*y(k,204)
         mat(k,1159) = .500_r8*rxt(k,321)*y(k,204)
         mat(k,1118) = .250_r8*rxt(k,450)*y(k,204)
         mat(k,1138) = .250_r8*rxt(k,455)*y(k,204)
         mat(k,1172) = .300_r8*rxt(k,391)*y(k,204)

         mat(k,447) = -(rxt(k,264)*y(k,223))
         mat(k,1611) = -rxt(k,264)*y(k,53)

         mat(k,1867) = rxt(k,261)*y(k,209)
         mat(k,1945) = rxt(k,261)*y(k,204)

         mat(k,1438) = -(rxt(k,176)*y(k,56) + rxt(k,232)*y(k,79) + rxt(k,265)*y(k,223) &
                      + (rxt(k,271) + rxt(k,272) + rxt(k,273)) * y(k,222))
         mat(k,2093) = -rxt(k,176)*y(k,54)
         mat(k,881) = -rxt(k,232)*y(k,54)
         mat(k,1695) = -rxt(k,265)*y(k,54)
         mat(k,1531) = -(rxt(k,271) + rxt(k,272) + rxt(k,273)) * y(k,54)

         mat(k,1015) = .100_r8*rxt(k,308)*y(k,140)
         mat(k,2155) = .100_r8*rxt(k,308)*y(k,29)

         mat(k,453) = -(rxt(k,228)*y(k,222) + rxt(k,245)*y(k,56) + rxt(k,246)*y(k,223))
         mat(k,1524) = -rxt(k,228)*y(k,55)
         mat(k,2076) = -rxt(k,245)*y(k,55)
         mat(k,1612) = -rxt(k,246)*y(k,55)

         mat(k,2105) = -(rxt(k,175)*y(k,42) + rxt(k,176)*y(k,54) + rxt(k,177)*y(k,83) &
                      + rxt(k,178)*y(k,85) + (rxt(k,179) + rxt(k,180)) * y(k,209) &
                      + rxt(k,181)*y(k,140) + rxt(k,188)*y(k,60) + rxt(k,197)*y(k,98) &
                      + rxt(k,238)*y(k,41) + rxt(k,240)*y(k,43) + rxt(k,243)*y(k,46) &
                      + rxt(k,245)*y(k,55) + rxt(k,286)*y(k,28))
         mat(k,2040) = -rxt(k,175)*y(k,56)
         mat(k,1447) = -rxt(k,176)*y(k,56)
         mat(k,1417) = -rxt(k,177)*y(k,56)
         mat(k,603) = -rxt(k,178)*y(k,56)
         mat(k,2016) = -(rxt(k,179) + rxt(k,180)) * y(k,56)
         mat(k,2167) = -rxt(k,181)*y(k,56)
         mat(k,938) = -rxt(k,188)*y(k,56)
         mat(k,811) = -rxt(k,197)*y(k,56)
         mat(k,477) = -rxt(k,238)*y(k,56)
         mat(k,597) = -rxt(k,240)*y(k,56)
         mat(k,376) = -rxt(k,243)*y(k,56)
         mat(k,457) = -rxt(k,245)*y(k,56)
         mat(k,297) = -rxt(k,286)*y(k,56)

         mat(k,1502) = rxt(k,216)*y(k,59)
         mat(k,107) = 4.000_r8*rxt(k,200)*y(k,222)
         mat(k,151) = rxt(k,201)*y(k,222)
         mat(k,122) = 2.000_r8*rxt(k,202)*y(k,222)
         mat(k,161) = 2.000_r8*rxt(k,203)*y(k,222)
         mat(k,126) = 2.000_r8*rxt(k,204)*y(k,222)
         mat(k,166) = rxt(k,205)*y(k,222)
         mat(k,130) = 2.000_r8*rxt(k,206)*y(k,222)
         mat(k,133) = 3.000_r8*rxt(k,242)*y(k,223)
         mat(k,376) = mat(k,376) + rxt(k,244)*y(k,223)
         mat(k,2066) = rxt(k,216)*y(k,19) + (4.000_r8*rxt(k,183)+2.000_r8*rxt(k,185)) &
                      *y(k,59) + rxt(k,187)*y(k,130) + rxt(k,192)*y(k,139)  &
                      + rxt(k,470)*y(k,156) + rxt(k,182)*y(k,204) + rxt(k,193) &
                      *y(k,223)
         mat(k,235) = rxt(k,237)*y(k,222)
         mat(k,231) = rxt(k,252)*y(k,222) + rxt(k,247)*y(k,223)
         mat(k,261) = rxt(k,253)*y(k,222) + rxt(k,248)*y(k,223)
         mat(k,314) = rxt(k,254)*y(k,222) + rxt(k,249)*y(k,223)
         mat(k,2190) = rxt(k,195)*y(k,139) + rxt(k,207)*y(k,222) + rxt(k,196)*y(k,223)
         mat(k,1858) = rxt(k,187)*y(k,59)
         mat(k,2265) = rxt(k,192)*y(k,59) + rxt(k,195)*y(k,91)
         mat(k,1245) = rxt(k,470)*y(k,59)
         mat(k,1909) = rxt(k,182)*y(k,59)
         mat(k,1543) = 4.000_r8*rxt(k,200)*y(k,33) + rxt(k,201)*y(k,34)  &
                      + 2.000_r8*rxt(k,202)*y(k,36) + 2.000_r8*rxt(k,203)*y(k,37)  &
                      + 2.000_r8*rxt(k,204)*y(k,38) + rxt(k,205)*y(k,39)  &
                      + 2.000_r8*rxt(k,206)*y(k,40) + rxt(k,237)*y(k,65) + rxt(k,252) &
                      *y(k,88) + rxt(k,253)*y(k,89) + rxt(k,254)*y(k,90) + rxt(k,207) &
                      *y(k,91)
         mat(k,1707) = 3.000_r8*rxt(k,242)*y(k,44) + rxt(k,244)*y(k,46) + rxt(k,193) &
                      *y(k,59) + rxt(k,247)*y(k,88) + rxt(k,248)*y(k,89) + rxt(k,249) &
                      *y(k,90) + rxt(k,196)*y(k,91)


         mat(k,2072) = rxt(k,188)*y(k,60)
         mat(k,2047) = 2.000_r8*rxt(k,184)*y(k,59)
         mat(k,931) = rxt(k,188)*y(k,56) + (rxt(k,535)+rxt(k,540)+rxt(k,545))*y(k,91)
         mat(k,2173) = (rxt(k,535)+rxt(k,540)+rxt(k,545))*y(k,60) + (rxt(k,530) &
                       +rxt(k,536)+rxt(k,541))*y(k,98)
         mat(k,807) = (rxt(k,530)+rxt(k,536)+rxt(k,541))*y(k,91)


         mat(k,2046) = 2.000_r8*rxt(k,209)*y(k,59)

         mat(k,2065) = -(rxt(k,182)*y(k,204) + (4._r8*rxt(k,183) + 4._r8*rxt(k,184) &
                      + 4._r8*rxt(k,185) + 4._r8*rxt(k,209)) * y(k,59) + rxt(k,186) &
                      *y(k,209) + rxt(k,187)*y(k,130) + rxt(k,189)*y(k,131) + rxt(k,192) &
                      *y(k,139) + (rxt(k,193) + rxt(k,194)) * y(k,223) + (rxt(k,215) &
                      + rxt(k,216) + rxt(k,217)) * y(k,19) + rxt(k,470)*y(k,156))
         mat(k,1908) = -rxt(k,182)*y(k,59)
         mat(k,2015) = -rxt(k,186)*y(k,59)
         mat(k,1857) = -rxt(k,187)*y(k,59)
         mat(k,2233) = -rxt(k,189)*y(k,59)
         mat(k,2264) = -rxt(k,192)*y(k,59)
         mat(k,1706) = -(rxt(k,193) + rxt(k,194)) * y(k,59)
         mat(k,1501) = -(rxt(k,215) + rxt(k,216) + rxt(k,217)) * y(k,59)
         mat(k,1244) = -rxt(k,470)*y(k,59)

         mat(k,2104) = rxt(k,197)*y(k,98) + rxt(k,181)*y(k,140) + rxt(k,180)*y(k,209)
         mat(k,937) = rxt(k,190)*y(k,139)
         mat(k,2189) = rxt(k,208)*y(k,222)
         mat(k,810) = rxt(k,197)*y(k,56) + rxt(k,198)*y(k,139) + rxt(k,199)*y(k,223)
         mat(k,2264) = mat(k,2264) + rxt(k,190)*y(k,60) + rxt(k,198)*y(k,98)
         mat(k,2166) = rxt(k,181)*y(k,56)
         mat(k,338) = rxt(k,475)*y(k,156)
         mat(k,1244) = mat(k,1244) + rxt(k,475)*y(k,142)
         mat(k,2015) = mat(k,2015) + rxt(k,180)*y(k,56)
         mat(k,1542) = rxt(k,208)*y(k,91)
         mat(k,1706) = mat(k,1706) + rxt(k,199)*y(k,98)

      end do

      end subroutine     nlnmat02

      subroutine     nlnmat03( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,933) = -(rxt(k,188)*y(k,56) + rxt(k,190)*y(k,139) + rxt(k,191)*y(k,223) &
                      + (rxt(k,535) + rxt(k,540) + rxt(k,545)) * y(k,91))
         mat(k,2084) = -rxt(k,188)*y(k,60)
         mat(k,2248) = -rxt(k,190)*y(k,60)
         mat(k,1664) = -rxt(k,191)*y(k,60)
         mat(k,2177) = -(rxt(k,535) + rxt(k,540) + rxt(k,545)) * y(k,60)

         mat(k,2052) = rxt(k,189)*y(k,131)
         mat(k,2213) = rxt(k,189)*y(k,59)


         mat(k,1088) = -(rxt(k,275)*y(k,223))
         mat(k,1676) = -rxt(k,275)*y(k,62)

         mat(k,981) = .230_r8*rxt(k,440)*y(k,140)
         mat(k,1422) = rxt(k,211)*y(k,42)
         mat(k,288) = .350_r8*rxt(k,277)*y(k,223)
         mat(k,557) = .630_r8*rxt(k,279)*y(k,140)
         mat(k,1011) = .560_r8*rxt(k,308)*y(k,140)
         mat(k,2025) = rxt(k,211)*y(k,17) + rxt(k,175)*y(k,56) + rxt(k,256)*y(k,132)  &
                      + rxt(k,257)*y(k,139) + rxt(k,258)*y(k,223)
         mat(k,372) = rxt(k,243)*y(k,56)
         mat(k,1225) = rxt(k,314)*y(k,132) + rxt(k,315)*y(k,223)
         mat(k,2089) = rxt(k,175)*y(k,42) + rxt(k,243)*y(k,46)
         mat(k,959) = rxt(k,302)*y(k,223)
         mat(k,856) = .620_r8*rxt(k,385)*y(k,140)
         mat(k,1213) = .650_r8*rxt(k,338)*y(k,140)
         mat(k,905) = .230_r8*rxt(k,443)*y(k,140)
         mat(k,1336) = .560_r8*rxt(k,352)*y(k,140)
         mat(k,1829) = .170_r8*rxt(k,411)*y(k,205) + .220_r8*rxt(k,336)*y(k,215)  &
                      + .400_r8*rxt(k,414)*y(k,216) + .350_r8*rxt(k,417)*y(k,218)  &
                      + .225_r8*rxt(k,452)*y(k,227) + .250_r8*rxt(k,393)*y(k,231)
         mat(k,1736) = rxt(k,256)*y(k,42) + rxt(k,314)*y(k,49) + .220_r8*rxt(k,335) &
                      *y(k,215) + .500_r8*rxt(k,394)*y(k,231)
         mat(k,2249) = rxt(k,257)*y(k,42) + rxt(k,464)*y(k,143)
         mat(k,2139) = .230_r8*rxt(k,440)*y(k,6) + .630_r8*rxt(k,279)*y(k,25)  &
                      + .560_r8*rxt(k,308)*y(k,29) + .620_r8*rxt(k,385)*y(k,104)  &
                      + .650_r8*rxt(k,338)*y(k,111) + .230_r8*rxt(k,443)*y(k,116)  &
                      + .560_r8*rxt(k,352)*y(k,117)
         mat(k,366) = rxt(k,464)*y(k,139) + rxt(k,465)*y(k,223)
         mat(k,1078) = .700_r8*rxt(k,461)*y(k,223)
         mat(k,1383) = .220_r8*rxt(k,332)*y(k,215) + .250_r8*rxt(k,390)*y(k,231)
         mat(k,1883) = .110_r8*rxt(k,333)*y(k,215) + .125_r8*rxt(k,450)*y(k,227)  &
                      + .200_r8*rxt(k,391)*y(k,231)
         mat(k,766) = .170_r8*rxt(k,411)*y(k,130) + .070_r8*rxt(k,410)*y(k,209)
         mat(k,1987) = .070_r8*rxt(k,410)*y(k,205) + .160_r8*rxt(k,413)*y(k,216)  &
                      + .140_r8*rxt(k,416)*y(k,218)
         mat(k,1314) = .220_r8*rxt(k,336)*y(k,130) + .220_r8*rxt(k,335)*y(k,132)  &
                      + .220_r8*rxt(k,332)*y(k,203) + .110_r8*rxt(k,333)*y(k,204)
         mat(k,729) = .400_r8*rxt(k,414)*y(k,130) + .160_r8*rxt(k,413)*y(k,209)
         mat(k,871) = .350_r8*rxt(k,417)*y(k,130) + .140_r8*rxt(k,416)*y(k,209)
         mat(k,1676) = mat(k,1676) + .350_r8*rxt(k,277)*y(k,24) + rxt(k,258)*y(k,42)  &
                      + rxt(k,315)*y(k,49) + rxt(k,302)*y(k,81) + rxt(k,465)*y(k,143)  &
                      + .700_r8*rxt(k,461)*y(k,184)
         mat(k,1123) = .225_r8*rxt(k,452)*y(k,130) + .125_r8*rxt(k,450)*y(k,204)
         mat(k,1176) = .250_r8*rxt(k,393)*y(k,130) + .500_r8*rxt(k,394)*y(k,132)  &
                      + .250_r8*rxt(k,390)*y(k,203) + .200_r8*rxt(k,391)*y(k,204)


         mat(k,971) = .270_r8*rxt(k,440)*y(k,140)
         mat(k,1006) = .200_r8*rxt(k,308)*y(k,140)
         mat(k,691) = rxt(k,295)*y(k,223)
         mat(k,617) = .500_r8*rxt(k,296)*y(k,223)
         mat(k,1087) = rxt(k,275)*y(k,223)
         mat(k,1051) = .800_r8*rxt(k,301)*y(k,223)
         mat(k,957) = rxt(k,302)*y(k,223)
         mat(k,943) = rxt(k,267)*y(k,223)
         mat(k,608) = .500_r8*rxt(k,351)*y(k,223)
         mat(k,896) = .270_r8*rxt(k,443)*y(k,140)
         mat(k,1332) = .100_r8*rxt(k,352)*y(k,140)
         mat(k,1816) = rxt(k,294)*y(k,203) + .900_r8*rxt(k,452)*y(k,227)
         mat(k,2124) = .270_r8*rxt(k,440)*y(k,6) + .200_r8*rxt(k,308)*y(k,29)  &
                      + .270_r8*rxt(k,443)*y(k,116) + .100_r8*rxt(k,352)*y(k,117)
         mat(k,1075) = 1.800_r8*rxt(k,461)*y(k,223)
         mat(k,1379) = rxt(k,294)*y(k,130) + 4.000_r8*rxt(k,291)*y(k,203)  &
                      + .900_r8*rxt(k,292)*y(k,204) + rxt(k,365)*y(k,211)  &
                      + 2.000_r8*rxt(k,341)*y(k,217) + rxt(k,390)*y(k,231)
         mat(k,1873) = .900_r8*rxt(k,292)*y(k,203) + rxt(k,342)*y(k,217)  &
                      + .500_r8*rxt(k,450)*y(k,227)
         mat(k,1975) = .450_r8*rxt(k,343)*y(k,217)
         mat(k,1253) = rxt(k,365)*y(k,203)
         mat(k,1358) = 2.000_r8*rxt(k,341)*y(k,203) + rxt(k,342)*y(k,204)  &
                      + .450_r8*rxt(k,343)*y(k,209) + 4.000_r8*rxt(k,344)*y(k,217)
         mat(k,1656) = rxt(k,295)*y(k,50) + .500_r8*rxt(k,296)*y(k,51) + rxt(k,275) &
                      *y(k,62) + .800_r8*rxt(k,301)*y(k,80) + rxt(k,302)*y(k,81)  &
                      + rxt(k,267)*y(k,93) + .500_r8*rxt(k,351)*y(k,115)  &
                      + 1.800_r8*rxt(k,461)*y(k,184)
         mat(k,1119) = .900_r8*rxt(k,452)*y(k,130) + .500_r8*rxt(k,450)*y(k,204)
         mat(k,1173) = rxt(k,390)*y(k,203)

         mat(k,253) = -(rxt(k,236)*y(k,222))
         mat(k,1521) = -rxt(k,236)*y(k,64)

         mat(k,148) = rxt(k,201)*y(k,222)
         mat(k,153) = rxt(k,227)*y(k,222)
         mat(k,159) = rxt(k,203)*y(k,222)
         mat(k,124) = 2.000_r8*rxt(k,204)*y(k,222)
         mat(k,163) = 2.000_r8*rxt(k,205)*y(k,222)
         mat(k,128) = rxt(k,206)*y(k,222)
         mat(k,112) = 2.000_r8*rxt(k,229)*y(k,222)
         mat(k,256) = rxt(k,253)*y(k,222) + rxt(k,248)*y(k,223)
         mat(k,309) = rxt(k,254)*y(k,222) + rxt(k,249)*y(k,223)
         mat(k,1521) = mat(k,1521) + rxt(k,201)*y(k,34) + rxt(k,227)*y(k,35)  &
                      + rxt(k,203)*y(k,37) + 2.000_r8*rxt(k,204)*y(k,38)  &
                      + 2.000_r8*rxt(k,205)*y(k,39) + rxt(k,206)*y(k,40)  &
                      + 2.000_r8*rxt(k,229)*y(k,84) + rxt(k,253)*y(k,89) + rxt(k,254) &
                      *y(k,90)
         mat(k,1584) = rxt(k,248)*y(k,89) + rxt(k,249)*y(k,90)

         mat(k,232) = -(rxt(k,237)*y(k,222))
         mat(k,1520) = -rxt(k,237)*y(k,65)

         mat(k,120) = rxt(k,202)*y(k,222)
         mat(k,158) = rxt(k,203)*y(k,222)
         mat(k,228) = rxt(k,252)*y(k,222) + rxt(k,247)*y(k,223)
         mat(k,1520) = mat(k,1520) + rxt(k,202)*y(k,36) + rxt(k,203)*y(k,37)  &
                      + rxt(k,252)*y(k,88)
         mat(k,1579) = rxt(k,247)*y(k,88)

         mat(k,200) = -(rxt(k,409)*y(k,223))
         mat(k,1573) = -rxt(k,409)*y(k,66)

         mat(k,194) = .180_r8*rxt(k,429)*y(k,223)
         mat(k,1573) = mat(k,1573) + .180_r8*rxt(k,429)*y(k,186)

         mat(k,303) = -(rxt(k,462)*y(k,132) + (rxt(k,463) + rxt(k,477)) * y(k,223))
         mat(k,1716) = -rxt(k,462)*y(k,67)
         mat(k,1591) = -(rxt(k,463) + rxt(k,477)) * y(k,67)












         mat(k,718) = rxt(k,297)*y(k,209)
         mat(k,1927) = rxt(k,297)*y(k,208)

         mat(k,879) = -(rxt(k,232)*y(k,54) + rxt(k,233)*y(k,83) + rxt(k,234)*y(k,235) &
                      + rxt(k,235)*y(k,95))
         mat(k,1435) = -rxt(k,232)*y(k,79)
         mat(k,1408) = -rxt(k,233)*y(k,79)
         mat(k,2275) = -rxt(k,234)*y(k,79)
         mat(k,1467) = -rxt(k,235)*y(k,79)

         mat(k,154) = rxt(k,227)*y(k,222)
         mat(k,164) = rxt(k,205)*y(k,222)
         mat(k,254) = 2.000_r8*rxt(k,236)*y(k,222)
         mat(k,233) = rxt(k,237)*y(k,222)
         mat(k,1528) = rxt(k,227)*y(k,35) + rxt(k,205)*y(k,39) + 2.000_r8*rxt(k,236) &
                      *y(k,64) + rxt(k,237)*y(k,65)

         mat(k,1053) = -(rxt(k,301)*y(k,223))
         mat(k,1673) = -rxt(k,301)*y(k,80)

         mat(k,579) = .700_r8*rxt(k,376)*y(k,223)
         mat(k,548) = .500_r8*rxt(k,377)*y(k,223)
         mat(k,381) = rxt(k,388)*y(k,223)
         mat(k,1826) = .050_r8*rxt(k,374)*y(k,212) + .530_r8*rxt(k,336)*y(k,215)  &
                      + .225_r8*rxt(k,452)*y(k,227) + .250_r8*rxt(k,393)*y(k,231)
         mat(k,1733) = .050_r8*rxt(k,375)*y(k,212) + .530_r8*rxt(k,335)*y(k,215)  &
                      + .250_r8*rxt(k,394)*y(k,231)
         mat(k,1382) = .530_r8*rxt(k,332)*y(k,215) + .250_r8*rxt(k,390)*y(k,231)
         mat(k,1880) = .260_r8*rxt(k,333)*y(k,215) + .125_r8*rxt(k,450)*y(k,227)  &
                      + .100_r8*rxt(k,391)*y(k,231)
         mat(k,1287) = .050_r8*rxt(k,374)*y(k,130) + .050_r8*rxt(k,375)*y(k,132)
         mat(k,1313) = .530_r8*rxt(k,336)*y(k,130) + .530_r8*rxt(k,335)*y(k,132)  &
                      + .530_r8*rxt(k,332)*y(k,203) + .260_r8*rxt(k,333)*y(k,204)
         mat(k,1673) = mat(k,1673) + .700_r8*rxt(k,376)*y(k,105) + .500_r8*rxt(k,377) &
                      *y(k,106) + rxt(k,388)*y(k,121)
         mat(k,1121) = .225_r8*rxt(k,452)*y(k,130) + .125_r8*rxt(k,450)*y(k,204)
         mat(k,1175) = .250_r8*rxt(k,393)*y(k,130) + .250_r8*rxt(k,394)*y(k,132)  &
                      + .250_r8*rxt(k,390)*y(k,203) + .100_r8*rxt(k,391)*y(k,204)

         mat(k,958) = -(rxt(k,302)*y(k,223))
         mat(k,1667) = -rxt(k,302)*y(k,81)

         mat(k,287) = .650_r8*rxt(k,277)*y(k,223)
         mat(k,1052) = .200_r8*rxt(k,301)*y(k,223)
         mat(k,1029) = rxt(k,389)*y(k,223)
         mat(k,1822) = rxt(k,400)*y(k,197) + .050_r8*rxt(k,374)*y(k,212)  &
                      + .400_r8*rxt(k,414)*y(k,216) + .170_r8*rxt(k,417)*y(k,218)  &
                      + .700_r8*rxt(k,420)*y(k,224) + .600_r8*rxt(k,427)*y(k,229)  &
                      + .250_r8*rxt(k,393)*y(k,231) + .340_r8*rxt(k,433)*y(k,232)  &
                      + .170_r8*rxt(k,436)*y(k,234)
         mat(k,1727) = .050_r8*rxt(k,375)*y(k,212) + .250_r8*rxt(k,394)*y(k,231)
         mat(k,499) = rxt(k,400)*y(k,130)
         mat(k,1380) = .250_r8*rxt(k,390)*y(k,231)
         mat(k,1876) = .100_r8*rxt(k,391)*y(k,231)
         mat(k,1981) = .160_r8*rxt(k,413)*y(k,216) + .070_r8*rxt(k,416)*y(k,218)
         mat(k,1286) = .050_r8*rxt(k,374)*y(k,130) + .050_r8*rxt(k,375)*y(k,132)
         mat(k,728) = .400_r8*rxt(k,414)*y(k,130) + .160_r8*rxt(k,413)*y(k,209)
         mat(k,870) = .170_r8*rxt(k,417)*y(k,130) + .070_r8*rxt(k,416)*y(k,209)
         mat(k,1667) = mat(k,1667) + .650_r8*rxt(k,277)*y(k,24) + .200_r8*rxt(k,301) &
                      *y(k,80) + rxt(k,389)*y(k,122)
         mat(k,461) = .700_r8*rxt(k,420)*y(k,130)
         mat(k,741) = .600_r8*rxt(k,427)*y(k,130)
         mat(k,1174) = .250_r8*rxt(k,393)*y(k,130) + .250_r8*rxt(k,394)*y(k,132)  &
                      + .250_r8*rxt(k,390)*y(k,203) + .100_r8*rxt(k,391)*y(k,204)
         mat(k,757) = .340_r8*rxt(k,433)*y(k,130)
         mat(k,506) = .170_r8*rxt(k,436)*y(k,130)

         mat(k,1453) = -((rxt(k,135) + rxt(k,136) + rxt(k,137)) * y(k,209) + rxt(k,141) &
                      *y(k,140))
         mat(k,2005) = -(rxt(k,135) + rxt(k,136) + rxt(k,137)) * y(k,82)
         mat(k,2156) = -rxt(k,141)*y(k,82)

         mat(k,2029) = rxt(k,258)*y(k,223)
         mat(k,1439) = rxt(k,272)*y(k,222)
         mat(k,2094) = rxt(k,177)*y(k,83)
         mat(k,882) = rxt(k,233)*y(k,83)
         mat(k,1411) = rxt(k,177)*y(k,56) + rxt(k,233)*y(k,79) + rxt(k,133)*y(k,139)  &
                      + rxt(k,125)*y(k,222) + rxt(k,142)*y(k,223)
         mat(k,817) = rxt(k,231)*y(k,222)
         mat(k,2179) = rxt(k,208)*y(k,222)
         mat(k,482) = rxt(k,163)*y(k,223)
         mat(k,2254) = rxt(k,133)*y(k,83) + rxt(k,145)*y(k,223)
         mat(k,368) = rxt(k,465)*y(k,223)
         mat(k,519) = rxt(k,471)*y(k,223)
         mat(k,1239) = rxt(k,476)*y(k,223)
         mat(k,1532) = rxt(k,272)*y(k,54) + rxt(k,125)*y(k,83) + rxt(k,231)*y(k,87)  &
                      + rxt(k,208)*y(k,91)
         mat(k,1696) = rxt(k,258)*y(k,42) + rxt(k,142)*y(k,83) + rxt(k,163)*y(k,118)  &
                      + rxt(k,145)*y(k,139) + rxt(k,465)*y(k,143) + rxt(k,471) &
                      *y(k,154) + rxt(k,476)*y(k,156)

         mat(k,1409) = -(rxt(k,125)*y(k,222) + rxt(k,133)*y(k,139) + rxt(k,142) &
                      *y(k,223) + rxt(k,177)*y(k,56) + rxt(k,233)*y(k,79))
         mat(k,1529) = -rxt(k,125)*y(k,83)
         mat(k,2251) = -rxt(k,133)*y(k,83)
         mat(k,1693) = -rxt(k,142)*y(k,83)
         mat(k,2091) = -rxt(k,177)*y(k,83)
         mat(k,880) = -rxt(k,233)*y(k,83)

         mat(k,1437) = rxt(k,273)*y(k,222)
         mat(k,1451) = rxt(k,135)*y(k,209)
         mat(k,2002) = rxt(k,135)*y(k,82)
         mat(k,1529) = mat(k,1529) + rxt(k,273)*y(k,54)

      end do

      end subroutine     nlnmat03

      subroutine     nlnmat04( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,111) = -(rxt(k,229)*y(k,222))
         mat(k,1509) = -rxt(k,229)*y(k,84)

         mat(k,600) = -(rxt(k,134)*y(k,139) + rxt(k,143)*y(k,223) + rxt(k,178)*y(k,56))
         mat(k,2243) = -rxt(k,134)*y(k,85)
         mat(k,1632) = -rxt(k,143)*y(k,85)
         mat(k,2080) = -rxt(k,178)*y(k,85)

         mat(k,1954) = 2.000_r8*rxt(k,149)*y(k,209)
         mat(k,1632) = mat(k,1632) + 2.000_r8*rxt(k,148)*y(k,223)


         mat(k,266) = rxt(k,478)*y(k,235)
         mat(k,2271) = rxt(k,478)*y(k,158)

         mat(k,815) = -(rxt(k,224)*y(k,139) + rxt(k,225)*y(k,223) + (rxt(k,230) &
                      + rxt(k,231)) * y(k,222))
         mat(k,2246) = -rxt(k,224)*y(k,87)
         mat(k,1654) = -rxt(k,225)*y(k,87)
         mat(k,1527) = -(rxt(k,230) + rxt(k,231)) * y(k,87)

         mat(k,1421) = rxt(k,211)*y(k,42) + rxt(k,212)*y(k,209)
         mat(k,2023) = rxt(k,211)*y(k,17)
         mat(k,1973) = rxt(k,212)*y(k,17)

         mat(k,227) = -(rxt(k,247)*y(k,223) + rxt(k,252)*y(k,222))
         mat(k,1578) = -rxt(k,247)*y(k,88)
         mat(k,1519) = -rxt(k,252)*y(k,88)

         mat(k,257) = -(rxt(k,248)*y(k,223) + rxt(k,253)*y(k,222))
         mat(k,1585) = -rxt(k,248)*y(k,89)
         mat(k,1522) = -rxt(k,253)*y(k,89)

         mat(k,310) = -(rxt(k,249)*y(k,223) + rxt(k,254)*y(k,222))
         mat(k,1592) = -rxt(k,249)*y(k,90)
         mat(k,1523) = -rxt(k,254)*y(k,90)

         mat(k,2192) = -(rxt(k,195)*y(k,139) + rxt(k,196)*y(k,223) + (rxt(k,207) &
                      + rxt(k,208)) * y(k,222) + (rxt(k,530) + rxt(k,536) + rxt(k,541) &
                      ) * y(k,98) + (rxt(k,535) + rxt(k,540) + rxt(k,545)) * y(k,60) &
                      + (rxt(k,537) + rxt(k,542)) * y(k,97))
         mat(k,2267) = -rxt(k,195)*y(k,91)
         mat(k,1709) = -rxt(k,196)*y(k,91)
         mat(k,1545) = -(rxt(k,207) + rxt(k,208)) * y(k,91)
         mat(k,812) = -(rxt(k,530) + rxt(k,536) + rxt(k,541)) * y(k,91)
         mat(k,939) = -(rxt(k,535) + rxt(k,540) + rxt(k,545)) * y(k,91)
         mat(k,789) = -(rxt(k,537) + rxt(k,542)) * y(k,91)

         mat(k,298) = rxt(k,286)*y(k,56)
         mat(k,478) = rxt(k,238)*y(k,56)
         mat(k,2042) = rxt(k,175)*y(k,56)
         mat(k,598) = rxt(k,240)*y(k,56)
         mat(k,377) = 2.000_r8*rxt(k,243)*y(k,56)
         mat(k,1448) = rxt(k,176)*y(k,56)
         mat(k,458) = rxt(k,245)*y(k,56)
         mat(k,2107) = rxt(k,286)*y(k,28) + rxt(k,238)*y(k,41) + rxt(k,175)*y(k,42)  &
                      + rxt(k,240)*y(k,43) + 2.000_r8*rxt(k,243)*y(k,46) + rxt(k,176) &
                      *y(k,54) + rxt(k,245)*y(k,55) + rxt(k,177)*y(k,83) + rxt(k,178) &
                      *y(k,85) + rxt(k,197)*y(k,98) + rxt(k,179)*y(k,209)
         mat(k,2068) = rxt(k,194)*y(k,223)
         mat(k,1418) = rxt(k,177)*y(k,56)
         mat(k,604) = rxt(k,178)*y(k,56)
         mat(k,812) = mat(k,812) + rxt(k,197)*y(k,56)
         mat(k,2018) = rxt(k,179)*y(k,56)
         mat(k,1709) = mat(k,1709) + rxt(k,194)*y(k,59)

         mat(k,185) = -(rxt(k,266)*y(k,223) + rxt(k,274)*y(k,222))
         mat(k,1571) = -rxt(k,266)*y(k,92)
         mat(k,1517) = -rxt(k,274)*y(k,92)

         mat(k,944) = -(rxt(k,267)*y(k,223))
         mat(k,1665) = -rxt(k,267)*y(k,93)

         mat(k,973) = .050_r8*rxt(k,440)*y(k,140)
         mat(k,286) = .350_r8*rxt(k,277)*y(k,223)
         mat(k,556) = .370_r8*rxt(k,279)*y(k,140)
         mat(k,1008) = .120_r8*rxt(k,308)*y(k,140)
         mat(k,854) = .110_r8*rxt(k,385)*y(k,140)
         mat(k,1212) = .330_r8*rxt(k,338)*y(k,140)
         mat(k,899) = .050_r8*rxt(k,443)*y(k,140)
         mat(k,1333) = .120_r8*rxt(k,352)*y(k,140)
         mat(k,1820) = rxt(k,270)*y(k,210)
         mat(k,2129) = .050_r8*rxt(k,440)*y(k,6) + .370_r8*rxt(k,279)*y(k,25)  &
                      + .120_r8*rxt(k,308)*y(k,29) + .110_r8*rxt(k,385)*y(k,104)  &
                      + .330_r8*rxt(k,338)*y(k,111) + .050_r8*rxt(k,443)*y(k,116)  &
                      + .120_r8*rxt(k,352)*y(k,117)
         mat(k,1979) = rxt(k,268)*y(k,210)
         mat(k,442) = rxt(k,270)*y(k,130) + rxt(k,268)*y(k,209)
         mat(k,1665) = mat(k,1665) + .350_r8*rxt(k,277)*y(k,24)


         mat(k,1433) = rxt(k,232)*y(k,79)
         mat(k,878) = rxt(k,232)*y(k,54) + rxt(k,233)*y(k,83) + rxt(k,235)*y(k,95)  &
                      + rxt(k,234)*y(k,235)
         mat(k,1407) = rxt(k,233)*y(k,79)
         mat(k,1466) = rxt(k,235)*y(k,79)
         mat(k,2273) = rxt(k,234)*y(k,79)

         mat(k,1471) = -(rxt(k,172)*y(k,223) + rxt(k,235)*y(k,79))
         mat(k,1697) = -rxt(k,172)*y(k,95)
         mat(k,883) = -rxt(k,235)*y(k,95)

         mat(k,2030) = rxt(k,256)*y(k,132)
         mat(k,1045) = rxt(k,288)*y(k,132)
         mat(k,1228) = rxt(k,314)*y(k,132)
         mat(k,934) = (rxt(k,535)+rxt(k,540)+rxt(k,545))*y(k,91)
         mat(k,305) = rxt(k,462)*y(k,132)
         mat(k,2180) = (rxt(k,535)+rxt(k,540)+rxt(k,545))*y(k,60)
         mat(k,2224) = rxt(k,171)*y(k,223)
         mat(k,1756) = rxt(k,256)*y(k,42) + rxt(k,288)*y(k,45) + rxt(k,314)*y(k,49)  &
                      + rxt(k,462)*y(k,67)
         mat(k,1697) = mat(k,1697) + rxt(k,171)*y(k,131)

         mat(k,421) = -(rxt(k,150)*y(k,223))
         mat(k,1608) = -rxt(k,150)*y(k,96)

         mat(k,2199) = rxt(k,169)*y(k,209)
         mat(k,1941) = rxt(k,169)*y(k,131)

         mat(k,784) = -(rxt(k,226)*y(k,139) + (rxt(k,537) + rxt(k,542)) * y(k,91))
         mat(k,2244) = -rxt(k,226)*y(k,97)
         mat(k,2175) = -(rxt(k,537) + rxt(k,542)) * y(k,97)

         mat(k,1487) = rxt(k,218)*y(k,209)
         mat(k,1970) = rxt(k,218)*y(k,19)

         mat(k,808) = -(rxt(k,197)*y(k,56) + rxt(k,198)*y(k,139) + rxt(k,199)*y(k,223) &
                      + (rxt(k,530) + rxt(k,536) + rxt(k,541)) * y(k,91))
         mat(k,2082) = -rxt(k,197)*y(k,98)
         mat(k,2245) = -rxt(k,198)*y(k,98)
         mat(k,1653) = -rxt(k,199)*y(k,98)
         mat(k,2176) = -(rxt(k,530) + rxt(k,536) + rxt(k,541)) * y(k,98)

         mat(k,2050) = rxt(k,186)*y(k,209)
         mat(k,932) = rxt(k,191)*y(k,223)
         mat(k,1972) = rxt(k,186)*y(k,59)
         mat(k,1653) = mat(k,1653) + rxt(k,191)*y(k,60)

         mat(k,1097) = -(rxt(k,331)*y(k,223))
         mat(k,1677) = -rxt(k,331)*y(k,99)

         mat(k,580) = .300_r8*rxt(k,376)*y(k,223)
         mat(k,549) = .500_r8*rxt(k,377)*y(k,223)
         mat(k,1830) = rxt(k,330)*y(k,206) + rxt(k,337)*y(k,215)
         mat(k,573) = rxt(k,330)*y(k,130)
         mat(k,1315) = rxt(k,337)*y(k,130)
         mat(k,1677) = mat(k,1677) + .300_r8*rxt(k,376)*y(k,105) + .500_r8*rxt(k,377) &
                      *y(k,106)

         mat(k,236) = -(rxt(k,362)*y(k,223))
         mat(k,1580) = -rxt(k,362)*y(k,100)

         mat(k,1110) = -(rxt(k,316)*y(k,223))
         mat(k,1678) = -rxt(k,316)*y(k,101)

         mat(k,581) = .700_r8*rxt(k,376)*y(k,223)
         mat(k,550) = .500_r8*rxt(k,377)*y(k,223)
         mat(k,609) = .500_r8*rxt(k,351)*y(k,223)
         mat(k,1831) = .050_r8*rxt(k,374)*y(k,212) + .220_r8*rxt(k,336)*y(k,215)  &
                      + .250_r8*rxt(k,393)*y(k,231)
         mat(k,1738) = .050_r8*rxt(k,375)*y(k,212) + .220_r8*rxt(k,335)*y(k,215)  &
                      + .250_r8*rxt(k,394)*y(k,231)
         mat(k,541) = .500_r8*rxt(k,320)*y(k,223)
         mat(k,1384) = .220_r8*rxt(k,332)*y(k,215) + .250_r8*rxt(k,390)*y(k,231)
         mat(k,1884) = .230_r8*rxt(k,333)*y(k,215) + .200_r8*rxt(k,321)*y(k,226)  &
                      + .100_r8*rxt(k,391)*y(k,231)
         mat(k,1290) = .050_r8*rxt(k,374)*y(k,130) + .050_r8*rxt(k,375)*y(k,132)
         mat(k,1316) = .220_r8*rxt(k,336)*y(k,130) + .220_r8*rxt(k,335)*y(k,132)  &
                      + .220_r8*rxt(k,332)*y(k,203) + .230_r8*rxt(k,333)*y(k,204)
         mat(k,1678) = mat(k,1678) + .700_r8*rxt(k,376)*y(k,105) + .500_r8*rxt(k,377) &
                      *y(k,106) + .500_r8*rxt(k,351)*y(k,115) + .500_r8*rxt(k,320) &
                      *y(k,152)
         mat(k,1160) = .200_r8*rxt(k,321)*y(k,204)
         mat(k,1177) = .250_r8*rxt(k,393)*y(k,130) + .250_r8*rxt(k,394)*y(k,132)  &
                      + .250_r8*rxt(k,390)*y(k,203) + .100_r8*rxt(k,391)*y(k,204)

         mat(k,332) = -(rxt(k,363)*y(k,223))
         mat(k,1596) = -rxt(k,363)*y(k,102)

         mat(k,1788) = .870_r8*rxt(k,374)*y(k,212)
         mat(k,1717) = .950_r8*rxt(k,375)*y(k,212)
         mat(k,1375) = rxt(k,370)*y(k,212)
         mat(k,1865) = .750_r8*rxt(k,371)*y(k,212)
         mat(k,1279) = .870_r8*rxt(k,374)*y(k,130) + .950_r8*rxt(k,375)*y(k,132)  &
                      + rxt(k,370)*y(k,203) + .750_r8*rxt(k,371)*y(k,204)

         mat(k,141) = -(rxt(k,364)*y(k,223))
         mat(k,1567) = -rxt(k,364)*y(k,103)

         mat(k,695) = .600_r8*rxt(k,387)*y(k,223)
         mat(k,1567) = mat(k,1567) + .600_r8*rxt(k,387)*y(k,109)

         mat(k,853) = -(rxt(k,378)*y(k,132) + rxt(k,385)*y(k,140) + rxt(k,386) &
                      *y(k,223))
         mat(k,1722) = -rxt(k,378)*y(k,104)
         mat(k,2126) = -rxt(k,385)*y(k,104)
         mat(k,1659) = -rxt(k,386)*y(k,104)

         mat(k,578) = -(rxt(k,376)*y(k,223))
         mat(k,1629) = -rxt(k,376)*y(k,105)

         mat(k,1802) = .080_r8*rxt(k,368)*y(k,211)
         mat(k,1250) = .080_r8*rxt(k,368)*y(k,130)

         mat(k,546) = -(rxt(k,377)*y(k,223))
         mat(k,1626) = -rxt(k,377)*y(k,106)

         mat(k,1800) = .080_r8*rxt(k,374)*y(k,212)
         mat(k,1280) = .080_r8*rxt(k,374)*y(k,130)

         mat(k,1198) = -(rxt(k,379)*y(k,203) + rxt(k,380)*y(k,204) + rxt(k,381) &
                      *y(k,209) + rxt(k,382)*y(k,130) + rxt(k,383)*y(k,132))
         mat(k,1386) = -rxt(k,379)*y(k,107)
         mat(k,1889) = -rxt(k,380)*y(k,107)
         mat(k,1993) = -rxt(k,381)*y(k,107)
         mat(k,1836) = -rxt(k,382)*y(k,107)
         mat(k,1743) = -rxt(k,383)*y(k,107)

         mat(k,857) = rxt(k,378)*y(k,132)
         mat(k,1743) = mat(k,1743) + rxt(k,378)*y(k,104)

         mat(k,391) = -(rxt(k,384)*y(k,223))
         mat(k,1604) = -rxt(k,384)*y(k,108)

         mat(k,1190) = rxt(k,381)*y(k,209)
         mat(k,1937) = rxt(k,381)*y(k,107)

         mat(k,696) = -(rxt(k,387)*y(k,223))
         mat(k,1643) = -rxt(k,387)*y(k,109)

         mat(k,1962) = rxt(k,367)*y(k,211) + rxt(k,372)*y(k,212)
         mat(k,1251) = rxt(k,367)*y(k,209)
         mat(k,1282) = rxt(k,372)*y(k,209)

         mat(k,80) = -(rxt(k,516)*y(k,223))
         mat(k,1559) = -rxt(k,516)*y(k,110)

         mat(k,1214) = -(rxt(k,338)*y(k,140) + rxt(k,339)*y(k,223))
         mat(k,2145) = -rxt(k,338)*y(k,111)
         mat(k,1684) = -rxt(k,339)*y(k,111)

         mat(k,858) = .300_r8*rxt(k,385)*y(k,140)
         mat(k,1837) = .360_r8*rxt(k,368)*y(k,211)
         mat(k,1744) = .400_r8*rxt(k,369)*y(k,211)
         mat(k,2145) = mat(k,2145) + .300_r8*rxt(k,385)*y(k,104)
         mat(k,1387) = .390_r8*rxt(k,365)*y(k,211)
         mat(k,1890) = .310_r8*rxt(k,366)*y(k,211)
         mat(k,1260) = .360_r8*rxt(k,368)*y(k,130) + .400_r8*rxt(k,369)*y(k,132)  &
                      + .390_r8*rxt(k,365)*y(k,203) + .310_r8*rxt(k,366)*y(k,204)

      end do

      end subroutine     nlnmat04

      subroutine     nlnmat05( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,316) = -(rxt(k,340)*y(k,223))
         mat(k,1593) = -rxt(k,340)*y(k,112)

         mat(k,1931) = rxt(k,334)*y(k,215)
         mat(k,1311) = rxt(k,334)*y(k,209)

         mat(k,512) = -(rxt(k,349)*y(k,223))
         mat(k,1621) = -rxt(k,349)*y(k,113)

         mat(k,1798) = .800_r8*rxt(k,358)*y(k,195)
         mat(k,916) = .800_r8*rxt(k,358)*y(k,130)

         mat(k,321) = -(rxt(k,350)*y(k,223))
         mat(k,1594) = -rxt(k,350)*y(k,114)

         mat(k,1932) = .800_r8*rxt(k,347)*y(k,219)
         mat(k,682) = .800_r8*rxt(k,347)*y(k,209)

         mat(k,607) = -(rxt(k,351)*y(k,223))
         mat(k,1633) = -rxt(k,351)*y(k,115)

         mat(k,2205) = rxt(k,354)*y(k,217)
         mat(k,1355) = rxt(k,354)*y(k,131)

         mat(k,897) = -(rxt(k,442)*y(k,132) + rxt(k,443)*y(k,140) + rxt(k,444) &
                      *y(k,223))
         mat(k,1723) = -rxt(k,442)*y(k,116)
         mat(k,2127) = -rxt(k,443)*y(k,116)
         mat(k,1662) = -rxt(k,444)*y(k,116)

         mat(k,1340) = -(rxt(k,352)*y(k,140) + rxt(k,353)*y(k,223))
         mat(k,2151) = -rxt(k,352)*y(k,117)
         mat(k,1690) = -rxt(k,353)*y(k,117)

         mat(k,861) = .200_r8*rxt(k,385)*y(k,140)
         mat(k,1842) = .560_r8*rxt(k,368)*y(k,211)
         mat(k,1750) = .600_r8*rxt(k,369)*y(k,211)
         mat(k,2151) = mat(k,2151) + .200_r8*rxt(k,385)*y(k,104)
         mat(k,1392) = .610_r8*rxt(k,365)*y(k,211)
         mat(k,1895) = .440_r8*rxt(k,366)*y(k,211)
         mat(k,1264) = .560_r8*rxt(k,368)*y(k,130) + .600_r8*rxt(k,369)*y(k,132)  &
                      + .610_r8*rxt(k,365)*y(k,203) + .440_r8*rxt(k,366)*y(k,204)

         mat(k,481) = -(rxt(k,151)*y(k,130) + (rxt(k,152) + rxt(k,153) + rxt(k,154) &
                      ) * y(k,131) + rxt(k,163)*y(k,223))
         mat(k,1795) = -rxt(k,151)*y(k,118)
         mat(k,2201) = -(rxt(k,152) + rxt(k,153) + rxt(k,154)) * y(k,118)
         mat(k,1617) = -rxt(k,163)*y(k,118)

         mat(k,189) = -((rxt(k,167) + rxt(k,168)) * y(k,222))
         mat(k,1518) = -(rxt(k,167) + rxt(k,168)) * y(k,119)

         mat(k,480) = rxt(k,152)*y(k,131)
         mat(k,2197) = rxt(k,152)*y(k,118)


         mat(k,2198) = rxt(k,170)*y(k,132)
         mat(k,1715) = rxt(k,170)*y(k,131)

         mat(k,379) = -(rxt(k,388)*y(k,223))
         mat(k,1602) = -rxt(k,388)*y(k,121)

         mat(k,1189) = .200_r8*rxt(k,380)*y(k,204)
         mat(k,1866) = .200_r8*rxt(k,380)*y(k,107)

         mat(k,1030) = -(rxt(k,389)*y(k,223))
         mat(k,1671) = -rxt(k,389)*y(k,122)

         mat(k,1194) = rxt(k,382)*y(k,130) + rxt(k,383)*y(k,132) + rxt(k,379)*y(k,203)  &
                      + .800_r8*rxt(k,380)*y(k,204)
         mat(k,1824) = rxt(k,382)*y(k,107)
         mat(k,1731) = rxt(k,383)*y(k,107)
         mat(k,1381) = rxt(k,379)*y(k,107)
         mat(k,1878) = .800_r8*rxt(k,380)*y(k,107)




         mat(k,102) = -(rxt(k,479)*y(k,223))
         mat(k,1563) = -rxt(k,479)*y(k,126)




         mat(k,1853) = -(rxt(k,151)*y(k,118) + rxt(k,160)*y(k,132) + rxt(k,164) &
                      *y(k,209) + rxt(k,165)*y(k,140) + rxt(k,166)*y(k,139) + rxt(k,187) &
                      *y(k,59) + rxt(k,219)*y(k,19) + rxt(k,262)*y(k,204) + rxt(k,270) &
                      *y(k,210) + rxt(k,283)*y(k,200) + rxt(k,294)*y(k,203) + rxt(k,298) &
                      *y(k,208) + rxt(k,311)*y(k,201) + rxt(k,319)*y(k,225) + rxt(k,323) &
                      *y(k,226) + (rxt(k,329) + rxt(k,330)) * y(k,206) + (rxt(k,336) &
                      + rxt(k,337)) * y(k,215) + rxt(k,345)*y(k,217) + rxt(k,348) &
                      *y(k,219) + (rxt(k,358) + rxt(k,359)) * y(k,195) + rxt(k,368) &
                      *y(k,211) + rxt(k,374)*y(k,212) + rxt(k,382)*y(k,107) + rxt(k,393) &
                      *y(k,231) + rxt(k,397)*y(k,194) + rxt(k,400)*y(k,197) + rxt(k,405) &
                      *y(k,199) + rxt(k,407)*y(k,202) + rxt(k,411)*y(k,205) + rxt(k,414) &
                      *y(k,216) + rxt(k,417)*y(k,218) + rxt(k,420)*y(k,224) + rxt(k,427) &
                      *y(k,229) + rxt(k,433)*y(k,232) + rxt(k,436)*y(k,234) + rxt(k,447) &
                      *y(k,221) + rxt(k,452)*y(k,227) + rxt(k,457)*y(k,228))
         mat(k,485) = -rxt(k,151)*y(k,130)
         mat(k,1761) = -rxt(k,160)*y(k,130)
         mat(k,2011) = -rxt(k,164)*y(k,130)
         mat(k,2162) = -rxt(k,165)*y(k,130)
         mat(k,2260) = -rxt(k,166)*y(k,130)
         mat(k,2061) = -rxt(k,187)*y(k,130)
         mat(k,1497) = -rxt(k,219)*y(k,130)
         mat(k,1904) = -rxt(k,262)*y(k,130)
         mat(k,443) = -rxt(k,270)*y(k,130)
         mat(k,829) = -rxt(k,283)*y(k,130)
         mat(k,1399) = -rxt(k,294)*y(k,130)
         mat(k,723) = -rxt(k,298)*y(k,130)
         mat(k,797) = -rxt(k,311)*y(k,130)
         mat(k,778) = -rxt(k,319)*y(k,130)
         mat(k,1165) = -rxt(k,323)*y(k,130)
         mat(k,574) = -(rxt(k,329) + rxt(k,330)) * y(k,130)
         mat(k,1325) = -(rxt(k,336) + rxt(k,337)) * y(k,130)
         mat(k,1367) = -rxt(k,345)*y(k,130)
         mat(k,687) = -rxt(k,348)*y(k,130)
         mat(k,926) = -(rxt(k,358) + rxt(k,359)) * y(k,130)
         mat(k,1270) = -rxt(k,368)*y(k,130)
         mat(k,1303) = -rxt(k,374)*y(k,130)
         mat(k,1206) = -rxt(k,382)*y(k,130)
         mat(k,1183) = -rxt(k,393)*y(k,130)
         mat(k,527) = -rxt(k,397)*y(k,130)
         mat(k,501) = -rxt(k,400)*y(k,130)
         mat(k,438) = -rxt(k,405)*y(k,130)
         mat(k,633) = -rxt(k,407)*y(k,130)
         mat(k,769) = -rxt(k,411)*y(k,130)
         mat(k,730) = -rxt(k,414)*y(k,130)
         mat(k,874) = -rxt(k,417)*y(k,130)
         mat(k,463) = -rxt(k,420)*y(k,130)
         mat(k,745) = -rxt(k,427)*y(k,130)
         mat(k,762) = -rxt(k,433)*y(k,130)
         mat(k,509) = -rxt(k,436)*y(k,130)
         mat(k,1070) = -rxt(k,447)*y(k,130)
         mat(k,1129) = -rxt(k,452)*y(k,130)
         mat(k,1151) = -rxt(k,457)*y(k,130)

         mat(k,485) = mat(k,485) + 2.000_r8*rxt(k,153)*y(k,131) + rxt(k,163)*y(k,223)
         mat(k,191) = 2.000_r8*rxt(k,167)*y(k,222)
         mat(k,2229) = 2.000_r8*rxt(k,153)*y(k,118) + rxt(k,156)*y(k,139) + rxt(k,472) &
                      *y(k,156)
         mat(k,2260) = mat(k,2260) + rxt(k,156)*y(k,131)
         mat(k,1242) = rxt(k,472)*y(k,131)
         mat(k,1538) = 2.000_r8*rxt(k,167)*y(k,119)
         mat(k,1702) = rxt(k,163)*y(k,118)

         mat(k,2237) = -((rxt(k,152) + rxt(k,153) + rxt(k,154)) * y(k,118) + (rxt(k,156) &
                      + rxt(k,158)) * y(k,139) + rxt(k,157)*y(k,140) + rxt(k,169) &
                      *y(k,209) + rxt(k,170)*y(k,132) + rxt(k,171)*y(k,223) + rxt(k,189) &
                      *y(k,59) + rxt(k,220)*y(k,19) + rxt(k,305)*y(k,203) + rxt(k,354) &
                      *y(k,217) + rxt(k,412)*y(k,205) + rxt(k,415)*y(k,216) + rxt(k,418) &
                      *y(k,218) + rxt(k,422)*y(k,147) + rxt(k,425)*y(k,194) + rxt(k,472) &
                      *y(k,156))
         mat(k,486) = -(rxt(k,152) + rxt(k,153) + rxt(k,154)) * y(k,131)
         mat(k,2268) = -(rxt(k,156) + rxt(k,158)) * y(k,131)
         mat(k,2170) = -rxt(k,157)*y(k,131)
         mat(k,2019) = -rxt(k,169)*y(k,131)
         mat(k,1769) = -rxt(k,170)*y(k,131)
         mat(k,1710) = -rxt(k,171)*y(k,131)
         mat(k,2069) = -rxt(k,189)*y(k,131)
         mat(k,1505) = -rxt(k,220)*y(k,131)
         mat(k,1404) = -rxt(k,305)*y(k,131)
         mat(k,1372) = -rxt(k,354)*y(k,131)
         mat(k,772) = -rxt(k,412)*y(k,131)
         mat(k,732) = -rxt(k,415)*y(k,131)
         mat(k,877) = -rxt(k,418)*y(k,131)
         mat(k,472) = -rxt(k,422)*y(k,131)
         mat(k,529) = -rxt(k,425)*y(k,131)
         mat(k,1247) = -rxt(k,472)*y(k,131)

         mat(k,647) = rxt(k,356)*y(k,223)
         mat(k,362) = rxt(k,327)*y(k,132)
         mat(k,1505) = mat(k,1505) + rxt(k,219)*y(k,130)
         mat(k,2069) = mat(k,2069) + rxt(k,187)*y(k,130)
         mat(k,425) = rxt(k,150)*y(k,223)
         mat(k,586) = .700_r8*rxt(k,376)*y(k,223)
         mat(k,1210) = rxt(k,382)*y(k,130) + rxt(k,383)*y(k,132)
         mat(k,1861) = rxt(k,219)*y(k,19) + rxt(k,187)*y(k,59) + rxt(k,382)*y(k,107)  &
                      + 2.000_r8*rxt(k,160)*y(k,132) + rxt(k,166)*y(k,139)  &
                      + rxt(k,165)*y(k,140) + rxt(k,397)*y(k,194) + rxt(k,358) &
                      *y(k,195) + rxt(k,400)*y(k,197) + rxt(k,405)*y(k,199)  &
                      + rxt(k,283)*y(k,200) + rxt(k,311)*y(k,201) + rxt(k,407) &
                      *y(k,202) + rxt(k,294)*y(k,203) + rxt(k,262)*y(k,204)  &
                      + rxt(k,411)*y(k,205) + rxt(k,329)*y(k,206) + rxt(k,298) &
                      *y(k,208) + rxt(k,164)*y(k,209) + rxt(k,270)*y(k,210)  &
                      + .920_r8*rxt(k,368)*y(k,211) + .920_r8*rxt(k,374)*y(k,212)  &
                      + rxt(k,336)*y(k,215) + rxt(k,414)*y(k,216) + rxt(k,345) &
                      *y(k,217) + rxt(k,417)*y(k,218) + rxt(k,348)*y(k,219)  &
                      + 1.600_r8*rxt(k,447)*y(k,221) + rxt(k,420)*y(k,224)  &
                      + rxt(k,319)*y(k,225) + rxt(k,323)*y(k,226) + .900_r8*rxt(k,452) &
                      *y(k,227) + .800_r8*rxt(k,457)*y(k,228) + rxt(k,427)*y(k,229)  &
                      + rxt(k,393)*y(k,231) + rxt(k,433)*y(k,232) + rxt(k,436) &
                      *y(k,234)
         mat(k,1769) = mat(k,1769) + rxt(k,327)*y(k,16) + rxt(k,383)*y(k,107)  &
                      + 2.000_r8*rxt(k,160)*y(k,130) + rxt(k,161)*y(k,139)  &
                      + rxt(k,159)*y(k,209) + rxt(k,369)*y(k,211) + rxt(k,375) &
                      *y(k,212) + rxt(k,335)*y(k,215) + rxt(k,346)*y(k,217)  &
                      + 2.000_r8*rxt(k,448)*y(k,221) + rxt(k,162)*y(k,223)  &
                      + rxt(k,394)*y(k,231)
         mat(k,842) = rxt(k,317)*y(k,223)
         mat(k,2268) = mat(k,2268) + rxt(k,166)*y(k,130) + rxt(k,161)*y(k,132)
         mat(k,2170) = mat(k,2170) + rxt(k,165)*y(k,130)
         mat(k,628) = rxt(k,454)*y(k,223)
         mat(k,529) = mat(k,529) + rxt(k,397)*y(k,130)
         mat(k,930) = rxt(k,358)*y(k,130)
         mat(k,503) = rxt(k,400)*y(k,130)
         mat(k,440) = rxt(k,405)*y(k,130)
         mat(k,833) = rxt(k,283)*y(k,130)
         mat(k,801) = rxt(k,311)*y(k,130)
         mat(k,636) = rxt(k,407)*y(k,130)
         mat(k,1404) = mat(k,1404) + rxt(k,294)*y(k,130)
         mat(k,1912) = rxt(k,262)*y(k,130) + .500_r8*rxt(k,445)*y(k,221)
         mat(k,772) = mat(k,772) + rxt(k,411)*y(k,130)
         mat(k,577) = rxt(k,329)*y(k,130)
         mat(k,726) = rxt(k,298)*y(k,130)
         mat(k,2019) = mat(k,2019) + rxt(k,164)*y(k,130) + rxt(k,159)*y(k,132)
         mat(k,446) = rxt(k,270)*y(k,130)
         mat(k,1275) = .920_r8*rxt(k,368)*y(k,130) + rxt(k,369)*y(k,132)
         mat(k,1308) = .920_r8*rxt(k,374)*y(k,130) + rxt(k,375)*y(k,132)
         mat(k,1329) = rxt(k,336)*y(k,130) + rxt(k,335)*y(k,132)
         mat(k,732) = mat(k,732) + rxt(k,414)*y(k,130)
         mat(k,1372) = mat(k,1372) + rxt(k,345)*y(k,130) + rxt(k,346)*y(k,132)
         mat(k,877) = mat(k,877) + rxt(k,417)*y(k,130)
         mat(k,689) = rxt(k,348)*y(k,130)
         mat(k,1074) = 1.600_r8*rxt(k,447)*y(k,130) + 2.000_r8*rxt(k,448)*y(k,132)  &
                      + .500_r8*rxt(k,445)*y(k,204)
         mat(k,1710) = mat(k,1710) + rxt(k,356)*y(k,1) + rxt(k,150)*y(k,96)  &
                      + .700_r8*rxt(k,376)*y(k,105) + rxt(k,162)*y(k,132) + rxt(k,317) &
                      *y(k,133) + rxt(k,454)*y(k,181)
         mat(k,465) = rxt(k,420)*y(k,130)
         mat(k,781) = rxt(k,319)*y(k,130)
         mat(k,1169) = rxt(k,323)*y(k,130)
         mat(k,1133) = .900_r8*rxt(k,452)*y(k,130)
         mat(k,1155) = .800_r8*rxt(k,457)*y(k,130)
         mat(k,747) = rxt(k,427)*y(k,130)
         mat(k,1187) = rxt(k,393)*y(k,130) + rxt(k,394)*y(k,132)
         mat(k,764) = rxt(k,433)*y(k,130)
         mat(k,511) = rxt(k,436)*y(k,130)

      end do

      end subroutine     nlnmat05

      subroutine     nlnmat06( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,1760) = -(rxt(k,159)*y(k,209) + rxt(k,160)*y(k,130) + rxt(k,161) &
                      *y(k,139) + rxt(k,162)*y(k,223) + rxt(k,170)*y(k,131) + rxt(k,256) &
                      *y(k,42) + rxt(k,288)*y(k,45) + rxt(k,307)*y(k,29) + rxt(k,314) &
                      *y(k,49) + rxt(k,327)*y(k,16) + rxt(k,335)*y(k,215) + rxt(k,346) &
                      *y(k,217) + rxt(k,369)*y(k,211) + rxt(k,375)*y(k,212) + rxt(k,378) &
                      *y(k,104) + rxt(k,383)*y(k,107) + rxt(k,394)*y(k,231) + rxt(k,439) &
                      *y(k,6) + rxt(k,442)*y(k,116) + rxt(k,448)*y(k,221) + rxt(k,459) &
                      *y(k,183) + rxt(k,462)*y(k,67))
         mat(k,2010) = -rxt(k,159)*y(k,132)
         mat(k,1852) = -rxt(k,160)*y(k,132)
         mat(k,2259) = -rxt(k,161)*y(k,132)
         mat(k,1701) = -rxt(k,162)*y(k,132)
         mat(k,2228) = -rxt(k,170)*y(k,132)
         mat(k,2034) = -rxt(k,256)*y(k,132)
         mat(k,1047) = -rxt(k,288)*y(k,132)
         mat(k,1017) = -rxt(k,307)*y(k,132)
         mat(k,1230) = -rxt(k,314)*y(k,132)
         mat(k,360) = -rxt(k,327)*y(k,132)
         mat(k,1324) = -rxt(k,335)*y(k,132)
         mat(k,1366) = -rxt(k,346)*y(k,132)
         mat(k,1269) = -rxt(k,369)*y(k,132)
         mat(k,1302) = -rxt(k,375)*y(k,132)
         mat(k,864) = -rxt(k,378)*y(k,132)
         mat(k,1205) = -rxt(k,383)*y(k,132)
         mat(k,1182) = -rxt(k,394)*y(k,132)
         mat(k,987) = -rxt(k,439)*y(k,132)
         mat(k,911) = -rxt(k,442)*y(k,132)
         mat(k,1069) = -rxt(k,448)*y(k,132)
         mat(k,1002) = -rxt(k,459)*y(k,132)
         mat(k,307) = -rxt(k,462)*y(k,132)

         mat(k,567) = rxt(k,221)*y(k,139)
         mat(k,2099) = rxt(k,188)*y(k,60)
         mat(k,936) = rxt(k,188)*y(k,56) + rxt(k,190)*y(k,139) + rxt(k,191)*y(k,223)
         mat(k,885) = rxt(k,235)*y(k,95)
         mat(k,1474) = rxt(k,235)*y(k,79) + rxt(k,172)*y(k,223)
         mat(k,612) = .500_r8*rxt(k,351)*y(k,223)
         mat(k,2228) = mat(k,2228) + rxt(k,158)*y(k,139) + rxt(k,157)*y(k,140)
         mat(k,2259) = mat(k,2259) + rxt(k,221)*y(k,20) + rxt(k,190)*y(k,60)  &
                      + rxt(k,158)*y(k,131)
         mat(k,2161) = rxt(k,157)*y(k,131)
         mat(k,534) = rxt(k,303)*y(k,223)
         mat(k,1701) = mat(k,1701) + rxt(k,191)*y(k,60) + rxt(k,172)*y(k,95)  &
                      + .500_r8*rxt(k,351)*y(k,115) + rxt(k,303)*y(k,145)

         mat(k,837) = -(rxt(k,317)*y(k,223))
         mat(k,1657) = -rxt(k,317)*y(k,133)

         mat(k,1007) = rxt(k,307)*y(k,132)
         mat(k,547) = .500_r8*rxt(k,377)*y(k,223)
         mat(k,393) = rxt(k,384)*y(k,223)
         mat(k,380) = rxt(k,388)*y(k,223)
         mat(k,1027) = rxt(k,389)*y(k,223)
         mat(k,1720) = rxt(k,307)*y(k,29)
         mat(k,1657) = mat(k,1657) + .500_r8*rxt(k,377)*y(k,106) + rxt(k,384)*y(k,108)  &
                      + rxt(k,388)*y(k,121) + rxt(k,389)*y(k,122)

         mat(k,385) = -(rxt(k,449)*y(k,223))
         mat(k,1603) = -rxt(k,449)*y(k,134)

         mat(k,1936) = rxt(k,446)*y(k,221)
         mat(k,1059) = rxt(k,446)*y(k,209)





         mat(k,2269) = -(rxt(k,130)*y(k,140) + 4._r8*rxt(k,131)*y(k,139) + rxt(k,133) &
                      *y(k,83) + rxt(k,134)*y(k,85) + rxt(k,139)*y(k,209) + rxt(k,145) &
                      *y(k,223) + (rxt(k,156) + rxt(k,158)) * y(k,131) + rxt(k,161) &
                      *y(k,132) + rxt(k,166)*y(k,130) + rxt(k,190)*y(k,60) + rxt(k,192) &
                      *y(k,59) + rxt(k,195)*y(k,91) + rxt(k,198)*y(k,98) + rxt(k,221) &
                      *y(k,20) + rxt(k,222)*y(k,19) + rxt(k,224)*y(k,87) + rxt(k,226) &
                      *y(k,97) + rxt(k,257)*y(k,42) + rxt(k,464)*y(k,143))
         mat(k,2171) = -rxt(k,130)*y(k,139)
         mat(k,1419) = -rxt(k,133)*y(k,139)
         mat(k,605) = -rxt(k,134)*y(k,139)
         mat(k,2020) = -rxt(k,139)*y(k,139)
         mat(k,1711) = -rxt(k,145)*y(k,139)
         mat(k,2238) = -(rxt(k,156) + rxt(k,158)) * y(k,139)
         mat(k,1770) = -rxt(k,161)*y(k,139)
         mat(k,1862) = -rxt(k,166)*y(k,139)
         mat(k,941) = -rxt(k,190)*y(k,139)
         mat(k,2070) = -rxt(k,192)*y(k,139)
         mat(k,2194) = -rxt(k,195)*y(k,139)
         mat(k,813) = -rxt(k,198)*y(k,139)
         mat(k,569) = -rxt(k,221)*y(k,139)
         mat(k,1506) = -rxt(k,222)*y(k,139)
         mat(k,821) = -rxt(k,224)*y(k,139)
         mat(k,790) = -rxt(k,226)*y(k,139)
         mat(k,2044) = -rxt(k,257)*y(k,139)
         mat(k,370) = -rxt(k,464)*y(k,139)

         mat(k,1464) = rxt(k,137)*y(k,209)
         mat(k,487) = rxt(k,151)*y(k,130) + rxt(k,152)*y(k,131)
         mat(k,1862) = mat(k,1862) + rxt(k,151)*y(k,118)
         mat(k,2238) = mat(k,2238) + rxt(k,152)*y(k,118)
         mat(k,2020) = mat(k,2020) + rxt(k,137)*y(k,82)
         mat(k,1711) = mat(k,1711) + 2.000_r8*rxt(k,147)*y(k,223)

         mat(k,2168) = -(rxt(k,129)*y(k,222) + rxt(k,130)*y(k,139) + rxt(k,140) &
                      *y(k,209) + rxt(k,141)*y(k,82) + rxt(k,146)*y(k,223) + rxt(k,157) &
                      *y(k,131) + rxt(k,165)*y(k,130) + rxt(k,181)*y(k,56) + rxt(k,213) &
                      *y(k,17) + rxt(k,279)*y(k,25) + rxt(k,308)*y(k,29) + rxt(k,338) &
                      *y(k,111) + rxt(k,352)*y(k,117) + rxt(k,385)*y(k,104) + rxt(k,423) &
                      *y(k,147) + rxt(k,440)*y(k,6) + rxt(k,443)*y(k,116) + rxt(k,468) &
                      *y(k,154) + rxt(k,474)*y(k,156))
         mat(k,1544) = -rxt(k,129)*y(k,140)
         mat(k,2266) = -rxt(k,130)*y(k,140)
         mat(k,2017) = -rxt(k,140)*y(k,140)
         mat(k,1462) = -rxt(k,141)*y(k,140)
         mat(k,1708) = -rxt(k,146)*y(k,140)
         mat(k,2235) = -rxt(k,157)*y(k,140)
         mat(k,1859) = -rxt(k,165)*y(k,140)
         mat(k,2106) = -rxt(k,181)*y(k,140)
         mat(k,1430) = -rxt(k,213)*y(k,140)
         mat(k,561) = -rxt(k,279)*y(k,140)
         mat(k,1022) = -rxt(k,308)*y(k,140)
         mat(k,1222) = -rxt(k,338)*y(k,140)
         mat(k,1351) = -rxt(k,352)*y(k,140)
         mat(k,868) = -rxt(k,385)*y(k,140)
         mat(k,471) = -rxt(k,423)*y(k,140)
         mat(k,992) = -rxt(k,440)*y(k,140)
         mat(k,914) = -rxt(k,443)*y(k,140)
         mat(k,521) = -rxt(k,468)*y(k,140)
         mat(k,1246) = -rxt(k,474)*y(k,140)

         mat(k,1403) = .150_r8*rxt(k,293)*y(k,209)
         mat(k,2017) = mat(k,2017) + .150_r8*rxt(k,293)*y(k,203) + .150_r8*rxt(k,343) &
                      *y(k,217)
         mat(k,1371) = .150_r8*rxt(k,343)*y(k,209)


         mat(k,335) = -(rxt(k,475)*y(k,156))
         mat(k,1234) = -rxt(k,475)*y(k,142)

         mat(k,1485) = rxt(k,215)*y(k,59)
         mat(k,2049) = rxt(k,215)*y(k,19) + 2.000_r8*rxt(k,185)*y(k,59)

         mat(k,363) = -(rxt(k,464)*y(k,139) + rxt(k,465)*y(k,223))
         mat(k,2240) = -rxt(k,464)*y(k,143)
         mat(k,1600) = -rxt(k,465)*y(k,143)


         mat(k,1092) = rxt(k,331)*y(k,223)
         mat(k,1785) = .100_r8*rxt(k,452)*y(k,227)
         mat(k,1582) = rxt(k,331)*y(k,99)
         mat(k,1116) = .100_r8*rxt(k,452)*y(k,130)

         mat(k,530) = -(rxt(k,303)*y(k,223))
         mat(k,1624) = -rxt(k,303)*y(k,145)

         mat(k,2203) = rxt(k,305)*y(k,203)
         mat(k,1376) = rxt(k,305)*y(k,131)


         mat(k,2196) = rxt(k,425)*y(k,194)
         mat(k,523) = rxt(k,425)*y(k,131)

         mat(k,469) = -(rxt(k,422)*y(k,131) + rxt(k,423)*y(k,140))
         mat(k,2200) = -rxt(k,422)*y(k,147)
         mat(k,2118) = -rxt(k,423)*y(k,147)

         mat(k,202) = .070_r8*rxt(k,409)*y(k,223)
         mat(k,1794) = rxt(k,407)*y(k,202)
         mat(k,180) = .060_r8*rxt(k,421)*y(k,223)
         mat(k,223) = .070_r8*rxt(k,437)*y(k,223)
         mat(k,630) = rxt(k,407)*y(k,130)
         mat(k,1615) = .070_r8*rxt(k,409)*y(k,66) + .060_r8*rxt(k,421)*y(k,148)  &
                      + .070_r8*rxt(k,437)*y(k,190)

         mat(k,178) = -(rxt(k,421)*y(k,223))
         mat(k,1570) = -rxt(k,421)*y(k,148)

         mat(k,170) = .530_r8*rxt(k,398)*y(k,223)
         mat(k,1570) = mat(k,1570) + .530_r8*rxt(k,398)*y(k,7)

         mat(k,340) = -(rxt(k,424)*y(k,223))
         mat(k,1597) = -rxt(k,424)*y(k,149)

         mat(k,1934) = rxt(k,419)*y(k,224)
         mat(k,459) = rxt(k,419)*y(k,209)



         mat(k,538) = -(rxt(k,320)*y(k,223))
         mat(k,1625) = -rxt(k,320)*y(k,152)

         mat(k,1953) = rxt(k,318)*y(k,225)
         mat(k,773) = rxt(k,318)*y(k,209)

         mat(k,403) = -(rxt(k,324)*y(k,223))
         mat(k,1606) = -rxt(k,324)*y(k,153)

         mat(k,1939) = .850_r8*rxt(k,322)*y(k,226)
         mat(k,1158) = .850_r8*rxt(k,322)*y(k,209)

         mat(k,517) = -(rxt(k,468)*y(k,140) + rxt(k,471)*y(k,223))
         mat(k,2119) = -rxt(k,468)*y(k,154)
         mat(k,1622) = -rxt(k,471)*y(k,154)


         mat(k,1237) = -(rxt(k,469)*y(k,19) + rxt(k,470)*y(k,59) + rxt(k,472)*y(k,131) &
                      + rxt(k,474)*y(k,140) + rxt(k,475)*y(k,142) + rxt(k,476) &
                      *y(k,223))
         mat(k,1489) = -rxt(k,469)*y(k,156)
         mat(k,2053) = -rxt(k,470)*y(k,156)
         mat(k,2218) = -rxt(k,472)*y(k,156)
         mat(k,2147) = -rxt(k,474)*y(k,156)
         mat(k,337) = -rxt(k,475)*y(k,156)
         mat(k,1686) = -rxt(k,476)*y(k,156)

         mat(k,2250) = rxt(k,464)*y(k,143)
         mat(k,2147) = mat(k,2147) + rxt(k,468)*y(k,154)
         mat(k,367) = rxt(k,464)*y(k,139)
         mat(k,518) = rxt(k,468)*y(k,140) + rxt(k,471)*y(k,223)
         mat(k,1686) = mat(k,1686) + rxt(k,471)*y(k,154)

         mat(k,844) = -(rxt(k,467)*y(k,223))
         mat(k,1658) = -rxt(k,467)*y(k,157)

         mat(k,1488) = rxt(k,469)*y(k,156)
         mat(k,2051) = rxt(k,470)*y(k,156)
         mat(k,304) = rxt(k,462)*y(k,132) + (rxt(k,463)+.500_r8*rxt(k,477))*y(k,223)
         mat(k,2211) = rxt(k,472)*y(k,156)
         mat(k,1721) = rxt(k,462)*y(k,67)
         mat(k,2125) = rxt(k,474)*y(k,156)
         mat(k,336) = rxt(k,475)*y(k,156)
         mat(k,365) = rxt(k,465)*y(k,223)
         mat(k,1236) = rxt(k,469)*y(k,19) + rxt(k,470)*y(k,59) + rxt(k,472)*y(k,131)  &
                      + rxt(k,474)*y(k,140) + rxt(k,475)*y(k,142) + rxt(k,476) &
                      *y(k,223)
         mat(k,1658) = mat(k,1658) + (rxt(k,463)+.500_r8*rxt(k,477))*y(k,67)  &
                      + rxt(k,465)*y(k,143) + rxt(k,476)*y(k,156)

         mat(k,267) = -(rxt(k,478)*y(k,235))
         mat(k,2272) = -rxt(k,478)*y(k,158)

         mat(k,843) = rxt(k,467)*y(k,223)
         mat(k,1587) = rxt(k,467)*y(k,157)















         mat(k,963) = .2202005_r8*rxt(k,504)*y(k,140)
         mat(k,888) = .0508005_r8*rxt(k,520)*y(k,140)
         mat(k,1772) = .1279005_r8*rxt(k,503)*y(k,196) + .0097005_r8*rxt(k,508) &
                      *y(k,198) + .0003005_r8*rxt(k,511)*y(k,213)  &
                      + .1056005_r8*rxt(k,515)*y(k,214) + .0245005_r8*rxt(k,519) &
                      *y(k,220) + .0154005_r8*rxt(k,525)*y(k,230)  &
                      + .0063005_r8*rxt(k,528)*y(k,233)
         mat(k,2111) = .2202005_r8*rxt(k,504)*y(k,6) + .0508005_r8*rxt(k,520)*y(k,116)
         mat(k,49) = .5931005_r8*rxt(k,522)*y(k,223)
         mat(k,55) = .1279005_r8*rxt(k,503)*y(k,130) + .2202005_r8*rxt(k,502)*y(k,209)
         mat(k,61) = .0097005_r8*rxt(k,508)*y(k,130) + .0023005_r8*rxt(k,507)*y(k,209)
         mat(k,1915) = .2202005_r8*rxt(k,502)*y(k,196) + .0023005_r8*rxt(k,507) &
                      *y(k,198) + .0031005_r8*rxt(k,510)*y(k,213)  &
                      + .2381005_r8*rxt(k,514)*y(k,214) + .0508005_r8*rxt(k,518) &
                      *y(k,220) + .1364005_r8*rxt(k,524)*y(k,230)  &
                      + .1677005_r8*rxt(k,527)*y(k,233)
         mat(k,67) = .0003005_r8*rxt(k,511)*y(k,130) + .0031005_r8*rxt(k,510)*y(k,209)
         mat(k,73) = .1056005_r8*rxt(k,515)*y(k,130) + .2381005_r8*rxt(k,514)*y(k,209)
         mat(k,81) = .0245005_r8*rxt(k,519)*y(k,130) + .0508005_r8*rxt(k,518)*y(k,209)
         mat(k,1549) = .5931005_r8*rxt(k,522)*y(k,178)
         mat(k,87) = .0154005_r8*rxt(k,525)*y(k,130) + .1364005_r8*rxt(k,524)*y(k,209)
         mat(k,93) = .0063005_r8*rxt(k,528)*y(k,130) + .1677005_r8*rxt(k,527)*y(k,209)

      end do

      end subroutine     nlnmat06

      subroutine     nlnmat07( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len

         mat(k,964) = .2067005_r8*rxt(k,504)*y(k,140)
         mat(k,889) = .1149005_r8*rxt(k,520)*y(k,140)
         mat(k,1773) = .1792005_r8*rxt(k,503)*y(k,196) + .0034005_r8*rxt(k,508) &
                      *y(k,198) + .0003005_r8*rxt(k,511)*y(k,213)  &
                      + .1026005_r8*rxt(k,515)*y(k,214) + .0082005_r8*rxt(k,519) &
                      *y(k,220) + .0452005_r8*rxt(k,525)*y(k,230)  &
                      + .0237005_r8*rxt(k,528)*y(k,233)
         mat(k,2112) = .2067005_r8*rxt(k,504)*y(k,6) + .1149005_r8*rxt(k,520)*y(k,116)
         mat(k,50) = .1534005_r8*rxt(k,522)*y(k,223)
         mat(k,56) = .1792005_r8*rxt(k,503)*y(k,130) + .2067005_r8*rxt(k,502)*y(k,209)
         mat(k,62) = .0034005_r8*rxt(k,508)*y(k,130) + .0008005_r8*rxt(k,507)*y(k,209)
         mat(k,1916) = .2067005_r8*rxt(k,502)*y(k,196) + .0008005_r8*rxt(k,507) &
                      *y(k,198) + .0035005_r8*rxt(k,510)*y(k,213)  &
                      + .1308005_r8*rxt(k,514)*y(k,214) + .1149005_r8*rxt(k,518) &
                      *y(k,220) + .0101005_r8*rxt(k,524)*y(k,230)  &
                      + .0174005_r8*rxt(k,527)*y(k,233)
         mat(k,68) = .0003005_r8*rxt(k,511)*y(k,130) + .0035005_r8*rxt(k,510)*y(k,209)
         mat(k,74) = .1026005_r8*rxt(k,515)*y(k,130) + .1308005_r8*rxt(k,514)*y(k,209)
         mat(k,82) = .0082005_r8*rxt(k,519)*y(k,130) + .1149005_r8*rxt(k,518)*y(k,209)
         mat(k,1550) = .1534005_r8*rxt(k,522)*y(k,178)
         mat(k,88) = .0452005_r8*rxt(k,525)*y(k,130) + .0101005_r8*rxt(k,524)*y(k,209)
         mat(k,94) = .0237005_r8*rxt(k,528)*y(k,130) + .0174005_r8*rxt(k,527)*y(k,209)


         mat(k,965) = .0653005_r8*rxt(k,504)*y(k,140)
         mat(k,890) = .0348005_r8*rxt(k,520)*y(k,140)
         mat(k,1774) = .0676005_r8*rxt(k,503)*y(k,196) + .1579005_r8*rxt(k,508) &
                      *y(k,198) + .0073005_r8*rxt(k,511)*y(k,213)  &
                      + .0521005_r8*rxt(k,515)*y(k,214) + .0772005_r8*rxt(k,519) &
                      *y(k,220) + .0966005_r8*rxt(k,525)*y(k,230)  &
                      + .0025005_r8*rxt(k,528)*y(k,233)
         mat(k,2113) = .0653005_r8*rxt(k,504)*y(k,6) + .0348005_r8*rxt(k,520)*y(k,116)
         mat(k,51) = .0459005_r8*rxt(k,522)*y(k,223)
         mat(k,57) = .0676005_r8*rxt(k,503)*y(k,130) + .0653005_r8*rxt(k,502)*y(k,209)
         mat(k,63) = .1579005_r8*rxt(k,508)*y(k,130) + .0843005_r8*rxt(k,507)*y(k,209)
         mat(k,1917) = .0653005_r8*rxt(k,502)*y(k,196) + .0843005_r8*rxt(k,507) &
                      *y(k,198) + .0003005_r8*rxt(k,510)*y(k,213)  &
                      + .0348005_r8*rxt(k,514)*y(k,214) + .0348005_r8*rxt(k,518) &
                      *y(k,220) + .0763005_r8*rxt(k,524)*y(k,230) + .086_r8*rxt(k,527) &
                      *y(k,233)
         mat(k,69) = .0073005_r8*rxt(k,511)*y(k,130) + .0003005_r8*rxt(k,510)*y(k,209)
         mat(k,75) = .0521005_r8*rxt(k,515)*y(k,130) + .0348005_r8*rxt(k,514)*y(k,209)
         mat(k,83) = .0772005_r8*rxt(k,519)*y(k,130) + .0348005_r8*rxt(k,518)*y(k,209)
         mat(k,1551) = .0459005_r8*rxt(k,522)*y(k,178)
         mat(k,89) = .0966005_r8*rxt(k,525)*y(k,130) + .0763005_r8*rxt(k,524)*y(k,209)
         mat(k,95) = .0025005_r8*rxt(k,528)*y(k,130) + .086_r8*rxt(k,527)*y(k,209)


         mat(k,966) = .1749305_r8*rxt(k,501)*y(k,132) + .1284005_r8*rxt(k,504) &
                      *y(k,140)
         mat(k,850) = .0590245_r8*rxt(k,509)*y(k,132) + .0033005_r8*rxt(k,512) &
                      *y(k,140)
         mat(k,891) = .1749305_r8*rxt(k,517)*y(k,132) + .0554005_r8*rxt(k,520) &
                      *y(k,140)
         mat(k,1775) = .079_r8*rxt(k,503)*y(k,196) + .0059005_r8*rxt(k,508)*y(k,198)  &
                      + .0057005_r8*rxt(k,511)*y(k,213) + .0143005_r8*rxt(k,515) &
                      *y(k,214) + .0332005_r8*rxt(k,519)*y(k,220)  &
                      + .0073005_r8*rxt(k,525)*y(k,230) + .011_r8*rxt(k,528)*y(k,233)
         mat(k,1713) = .1749305_r8*rxt(k,501)*y(k,6) + .0590245_r8*rxt(k,509)*y(k,104)  &
                      + .1749305_r8*rxt(k,517)*y(k,116)
         mat(k,2114) = .1284005_r8*rxt(k,504)*y(k,6) + .0033005_r8*rxt(k,512)*y(k,104)  &
                      + .0554005_r8*rxt(k,520)*y(k,116)
         mat(k,52) = .0085005_r8*rxt(k,522)*y(k,223)
         mat(k,58) = .079_r8*rxt(k,503)*y(k,130) + .1284005_r8*rxt(k,502)*y(k,209)
         mat(k,64) = .0059005_r8*rxt(k,508)*y(k,130) + .0443005_r8*rxt(k,507)*y(k,209)
         mat(k,1918) = .1284005_r8*rxt(k,502)*y(k,196) + .0443005_r8*rxt(k,507) &
                      *y(k,198) + .0271005_r8*rxt(k,510)*y(k,213)  &
                      + .0076005_r8*rxt(k,514)*y(k,214) + .0554005_r8*rxt(k,518) &
                      *y(k,220) + .2157005_r8*rxt(k,524)*y(k,230)  &
                      + .0512005_r8*rxt(k,527)*y(k,233)
         mat(k,70) = .0057005_r8*rxt(k,511)*y(k,130) + .0271005_r8*rxt(k,510)*y(k,209)
         mat(k,76) = .0143005_r8*rxt(k,515)*y(k,130) + .0076005_r8*rxt(k,514)*y(k,209)
         mat(k,84) = .0332005_r8*rxt(k,519)*y(k,130) + .0554005_r8*rxt(k,518)*y(k,209)
         mat(k,1552) = .0085005_r8*rxt(k,522)*y(k,178)
         mat(k,90) = .0073005_r8*rxt(k,525)*y(k,130) + .2157005_r8*rxt(k,524)*y(k,209)
         mat(k,96) = .011_r8*rxt(k,528)*y(k,130) + .0512005_r8*rxt(k,527)*y(k,209)


         mat(k,967) = .5901905_r8*rxt(k,501)*y(k,132) + .114_r8*rxt(k,504)*y(k,140)
         mat(k,851) = .0250245_r8*rxt(k,509)*y(k,132)
         mat(k,892) = .5901905_r8*rxt(k,517)*y(k,132) + .1278005_r8*rxt(k,520) &
                      *y(k,140)
         mat(k,1776) = .1254005_r8*rxt(k,503)*y(k,196) + .0536005_r8*rxt(k,508) &
                      *y(k,198) + .0623005_r8*rxt(k,511)*y(k,213)  &
                      + .0166005_r8*rxt(k,515)*y(k,214) + .130_r8*rxt(k,519)*y(k,220)  &
                      + .238_r8*rxt(k,525)*y(k,230) + .1185005_r8*rxt(k,528)*y(k,233)
         mat(k,1714) = .5901905_r8*rxt(k,501)*y(k,6) + .0250245_r8*rxt(k,509)*y(k,104)  &
                      + .5901905_r8*rxt(k,517)*y(k,116)
         mat(k,2115) = .114_r8*rxt(k,504)*y(k,6) + .1278005_r8*rxt(k,520)*y(k,116)
         mat(k,53) = .0128005_r8*rxt(k,522)*y(k,223)
         mat(k,59) = .1254005_r8*rxt(k,503)*y(k,130) + .114_r8*rxt(k,502)*y(k,209)
         mat(k,65) = .0536005_r8*rxt(k,508)*y(k,130) + .1621005_r8*rxt(k,507)*y(k,209)
         mat(k,1919) = .114_r8*rxt(k,502)*y(k,196) + .1621005_r8*rxt(k,507)*y(k,198)  &
                      + .0474005_r8*rxt(k,510)*y(k,213) + .0113005_r8*rxt(k,514) &
                      *y(k,214) + .1278005_r8*rxt(k,518)*y(k,220)  &
                      + .0738005_r8*rxt(k,524)*y(k,230) + .1598005_r8*rxt(k,527) &
                      *y(k,233)
         mat(k,71) = .0623005_r8*rxt(k,511)*y(k,130) + .0474005_r8*rxt(k,510)*y(k,209)
         mat(k,77) = .0166005_r8*rxt(k,515)*y(k,130) + .0113005_r8*rxt(k,514)*y(k,209)
         mat(k,85) = .130_r8*rxt(k,519)*y(k,130) + .1278005_r8*rxt(k,518)*y(k,209)
         mat(k,1553) = .0128005_r8*rxt(k,522)*y(k,178)
         mat(k,91) = .238_r8*rxt(k,525)*y(k,130) + .0738005_r8*rxt(k,524)*y(k,209)
         mat(k,97) = .1185005_r8*rxt(k,528)*y(k,130) + .1598005_r8*rxt(k,527)*y(k,209)


         mat(k,54) = -(rxt(k,522)*y(k,223))
         mat(k,1554) = -rxt(k,522)*y(k,178)


         mat(k,195) = .100_r8*rxt(k,429)*y(k,223)
         mat(k,213) = .230_r8*rxt(k,431)*y(k,223)
         mat(k,1574) = .100_r8*rxt(k,429)*y(k,186) + .230_r8*rxt(k,431)*y(k,188)

         mat(k,648) = -(rxt(k,453)*y(k,223))
         mat(k,1638) = -rxt(k,453)*y(k,180)

         mat(k,1957) = rxt(k,451)*y(k,227)
         mat(k,1117) = rxt(k,451)*y(k,209)

         mat(k,623) = -(rxt(k,454)*y(k,223))
         mat(k,1635) = -rxt(k,454)*y(k,181)

         mat(k,1804) = .200_r8*rxt(k,447)*y(k,221) + .200_r8*rxt(k,457)*y(k,228)
         mat(k,1868) = .500_r8*rxt(k,445)*y(k,221)
         mat(k,1060) = .200_r8*rxt(k,447)*y(k,130) + .500_r8*rxt(k,445)*y(k,204)
         mat(k,1137) = .200_r8*rxt(k,457)*y(k,130)

         mat(k,488) = -(rxt(k,458)*y(k,223))
         mat(k,1618) = -rxt(k,458)*y(k,182)

         mat(k,1949) = rxt(k,456)*y(k,228)
         mat(k,1136) = rxt(k,456)*y(k,209)

         mat(k,996) = -(rxt(k,459)*y(k,132) + rxt(k,460)*y(k,223))
         mat(k,1729) = -rxt(k,459)*y(k,183)
         mat(k,1669) = -rxt(k,460)*y(k,183)

         mat(k,977) = .330_r8*rxt(k,440)*y(k,140)
         mat(k,902) = .330_r8*rxt(k,443)*y(k,140)
         mat(k,1823) = .800_r8*rxt(k,447)*y(k,221) + .800_r8*rxt(k,457)*y(k,228)
         mat(k,1729) = mat(k,1729) + rxt(k,448)*y(k,221)
         mat(k,2133) = .330_r8*rxt(k,440)*y(k,6) + .330_r8*rxt(k,443)*y(k,116)
         mat(k,624) = rxt(k,454)*y(k,223)
         mat(k,1877) = .500_r8*rxt(k,445)*y(k,221) + rxt(k,455)*y(k,228)
         mat(k,1062) = .800_r8*rxt(k,447)*y(k,130) + rxt(k,448)*y(k,132)  &
                      + .500_r8*rxt(k,445)*y(k,204)
         mat(k,1669) = mat(k,1669) + rxt(k,454)*y(k,181)
         mat(k,1140) = .800_r8*rxt(k,457)*y(k,130) + rxt(k,455)*y(k,204)

         mat(k,1077) = -(rxt(k,461)*y(k,223))
         mat(k,1675) = -rxt(k,461)*y(k,184)

         mat(k,980) = .300_r8*rxt(k,440)*y(k,140)
         mat(k,904) = .300_r8*rxt(k,443)*y(k,140)
         mat(k,1828) = .900_r8*rxt(k,452)*y(k,227)
         mat(k,2138) = .300_r8*rxt(k,440)*y(k,6) + .300_r8*rxt(k,443)*y(k,116)
         mat(k,1882) = rxt(k,450)*y(k,227)
         mat(k,1122) = .900_r8*rxt(k,452)*y(k,130) + rxt(k,450)*y(k,204)

         mat(k,661) = -(rxt(k,428)*y(k,223))
         mat(k,1639) = -rxt(k,428)*y(k,185)

         mat(k,1958) = rxt(k,426)*y(k,229)
         mat(k,736) = rxt(k,426)*y(k,209)

         mat(k,193) = -(rxt(k,429)*y(k,223))
         mat(k,1572) = -rxt(k,429)*y(k,186)

         mat(k,209) = -(rxt(k,395)*y(k,223))
         mat(k,1575) = -rxt(k,395)*y(k,187)

         mat(k,1928) = rxt(k,392)*y(k,231)
         mat(k,1171) = rxt(k,392)*y(k,209)

         mat(k,214) = -(rxt(k,431)*y(k,223))
         mat(k,1576) = -rxt(k,431)*y(k,188)

         mat(k,707) = -(rxt(k,434)*y(k,223))
         mat(k,1644) = -rxt(k,434)*y(k,189)

         mat(k,1963) = rxt(k,432)*y(k,232)
         mat(k,752) = rxt(k,432)*y(k,209)

         mat(k,222) = -(rxt(k,437)*y(k,223))
         mat(k,1577) = -rxt(k,437)*y(k,190)

         mat(k,215) = .150_r8*rxt(k,431)*y(k,223)
         mat(k,1577) = mat(k,1577) + .150_r8*rxt(k,431)*y(k,188)

         mat(k,427) = -(rxt(k,438)*y(k,223))
         mat(k,1609) = -rxt(k,438)*y(k,191)

         mat(k,1942) = rxt(k,435)*y(k,234)
         mat(k,504) = rxt(k,435)*y(k,209)

         mat(k,524) = -(rxt(k,396)*y(k,209) + rxt(k,397)*y(k,130) + rxt(k,425) &
                      *y(k,131))
         mat(k,1952) = -rxt(k,396)*y(k,194)
         mat(k,1799) = -rxt(k,397)*y(k,194)
         mat(k,2202) = -rxt(k,425)*y(k,194)

         mat(k,245) = rxt(k,402)*y(k,223)
         mat(k,1623) = rxt(k,402)*y(k,22)

         mat(k,921) = -(rxt(k,357)*y(k,209) + (rxt(k,358) + rxt(k,359)) * y(k,130))
         mat(k,1978) = -rxt(k,357)*y(k,195)
         mat(k,1819) = -(rxt(k,358) + rxt(k,359)) * y(k,195)

         mat(k,675) = rxt(k,360)*y(k,223)
         mat(k,251) = rxt(k,361)*y(k,223)
         mat(k,1663) = rxt(k,360)*y(k,2) + rxt(k,361)*y(k,15)

         mat(k,60) = -(rxt(k,502)*y(k,209) + rxt(k,503)*y(k,130))
         mat(k,1920) = -rxt(k,502)*y(k,196)
         mat(k,1777) = -rxt(k,503)*y(k,196)

         mat(k,968) = rxt(k,505)*y(k,223)
         mat(k,1555) = rxt(k,505)*y(k,6)

         mat(k,497) = -(rxt(k,399)*y(k,209) + rxt(k,400)*y(k,130))
         mat(k,1950) = -rxt(k,399)*y(k,197)
         mat(k,1796) = -rxt(k,400)*y(k,197)

         mat(k,171) = .350_r8*rxt(k,398)*y(k,223)
         mat(k,399) = rxt(k,401)*y(k,223)
         mat(k,1619) = .350_r8*rxt(k,398)*y(k,7) + rxt(k,401)*y(k,8)

         mat(k,66) = -(rxt(k,507)*y(k,209) + rxt(k,508)*y(k,130))
         mat(k,1921) = -rxt(k,507)*y(k,198)
         mat(k,1778) = -rxt(k,508)*y(k,198)

         mat(k,167) = rxt(k,506)*y(k,223)
         mat(k,1556) = rxt(k,506)*y(k,7)

         mat(k,435) = -(rxt(k,403)*y(k,209) + rxt(k,405)*y(k,130))
         mat(k,1943) = -rxt(k,403)*y(k,199)
         mat(k,1790) = -rxt(k,405)*y(k,199)

         mat(k,328) = rxt(k,404)*y(k,223)
         mat(k,196) = .070_r8*rxt(k,429)*y(k,223)
         mat(k,216) = .060_r8*rxt(k,431)*y(k,223)
         mat(k,1610) = rxt(k,404)*y(k,23) + .070_r8*rxt(k,429)*y(k,186)  &
                      + .060_r8*rxt(k,431)*y(k,188)

         mat(k,826) = -(4._r8*rxt(k,280)*y(k,200) + rxt(k,281)*y(k,204) + rxt(k,282) &
                      *y(k,209) + rxt(k,283)*y(k,130))
         mat(k,1872) = -rxt(k,281)*y(k,200)
         mat(k,1974) = -rxt(k,282)*y(k,200)
         mat(k,1815) = -rxt(k,283)*y(k,200)

         mat(k,346) = .500_r8*rxt(k,285)*y(k,223)
         mat(k,295) = rxt(k,286)*y(k,56) + rxt(k,287)*y(k,223)
         mat(k,2083) = rxt(k,286)*y(k,28)
         mat(k,1655) = .500_r8*rxt(k,285)*y(k,27) + rxt(k,287)*y(k,28)

      end do

      end subroutine     nlnmat07

      subroutine     nlnmat08( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,793) = -(rxt(k,309)*y(k,204) + rxt(k,310)*y(k,209) + rxt(k,311) &
                      *y(k,130))
         mat(k,1870) = -rxt(k,309)*y(k,201)
         mat(k,1971) = -rxt(k,310)*y(k,201)
         mat(k,1814) = -rxt(k,311)*y(k,201)

         mat(k,410) = rxt(k,312)*y(k,223)
         mat(k,116) = rxt(k,313)*y(k,223)
         mat(k,1651) = rxt(k,312)*y(k,30) + rxt(k,313)*y(k,31)

         mat(k,631) = -(rxt(k,406)*y(k,209) + rxt(k,407)*y(k,130))
         mat(k,1956) = -rxt(k,406)*y(k,202)
         mat(k,1805) = -rxt(k,407)*y(k,202)

         mat(k,273) = rxt(k,408)*y(k,223)
         mat(k,1805) = mat(k,1805) + rxt(k,397)*y(k,194)
         mat(k,2121) = rxt(k,423)*y(k,147)
         mat(k,470) = rxt(k,423)*y(k,140)
         mat(k,525) = rxt(k,397)*y(k,130) + .400_r8*rxt(k,396)*y(k,209)
         mat(k,1956) = mat(k,1956) + .400_r8*rxt(k,396)*y(k,194)
         mat(k,1636) = rxt(k,408)*y(k,32)

         mat(k,1394) = -(4._r8*rxt(k,291)*y(k,203) + rxt(k,292)*y(k,204) + rxt(k,293) &
                      *y(k,209) + rxt(k,294)*y(k,130) + rxt(k,305)*y(k,131) + rxt(k,332) &
                      *y(k,215) + rxt(k,365)*y(k,211) + rxt(k,370)*y(k,212) + rxt(k,379) &
                      *y(k,107) + rxt(k,390)*y(k,231))
         mat(k,1897) = -rxt(k,292)*y(k,203)
         mat(k,2001) = -rxt(k,293)*y(k,203)
         mat(k,1844) = -rxt(k,294)*y(k,203)
         mat(k,2220) = -rxt(k,305)*y(k,203)
         mat(k,1321) = -rxt(k,332)*y(k,203)
         mat(k,1266) = -rxt(k,365)*y(k,203)
         mat(k,1299) = -rxt(k,370)*y(k,203)
         mat(k,1202) = -rxt(k,379)*y(k,203)
         mat(k,1180) = -rxt(k,390)*y(k,203)

         mat(k,985) = .060_r8*rxt(k,440)*y(k,140)
         mat(k,1044) = rxt(k,288)*y(k,132) + rxt(k,289)*y(k,223)
         mat(k,1227) = rxt(k,314)*y(k,132) + rxt(k,315)*y(k,223)
         mat(k,618) = .500_r8*rxt(k,296)*y(k,223)
         mat(k,862) = .080_r8*rxt(k,385)*y(k,140)
         mat(k,1218) = .100_r8*rxt(k,338)*y(k,140)
         mat(k,909) = .060_r8*rxt(k,443)*y(k,140)
         mat(k,1342) = .280_r8*rxt(k,352)*y(k,140)
         mat(k,1844) = mat(k,1844) + .530_r8*rxt(k,336)*y(k,215) + rxt(k,345)*y(k,217)  &
                      + rxt(k,348)*y(k,219) + rxt(k,323)*y(k,226)
         mat(k,1752) = rxt(k,288)*y(k,45) + rxt(k,314)*y(k,49) + .530_r8*rxt(k,335) &
                      *y(k,215) + rxt(k,346)*y(k,217)
         mat(k,2153) = .060_r8*rxt(k,440)*y(k,6) + .080_r8*rxt(k,385)*y(k,104)  &
                      + .100_r8*rxt(k,338)*y(k,111) + .060_r8*rxt(k,443)*y(k,116)  &
                      + .280_r8*rxt(k,352)*y(k,117)
         mat(k,1080) = .650_r8*rxt(k,461)*y(k,223)
         mat(k,1394) = mat(k,1394) + .530_r8*rxt(k,332)*y(k,215)
         mat(k,1897) = mat(k,1897) + .260_r8*rxt(k,333)*y(k,215) + rxt(k,342)*y(k,217)  &
                      + .300_r8*rxt(k,321)*y(k,226)
         mat(k,2001) = mat(k,2001) + .450_r8*rxt(k,343)*y(k,217) + .200_r8*rxt(k,347) &
                      *y(k,219) + .150_r8*rxt(k,322)*y(k,226)
         mat(k,1321) = mat(k,1321) + .530_r8*rxt(k,336)*y(k,130) + .530_r8*rxt(k,335) &
                      *y(k,132) + .530_r8*rxt(k,332)*y(k,203) + .260_r8*rxt(k,333) &
                      *y(k,204)
         mat(k,1363) = rxt(k,345)*y(k,130) + rxt(k,346)*y(k,132) + rxt(k,342)*y(k,204)  &
                      + .450_r8*rxt(k,343)*y(k,209) + 4.000_r8*rxt(k,344)*y(k,217)
         mat(k,685) = rxt(k,348)*y(k,130) + .200_r8*rxt(k,347)*y(k,209)
         mat(k,1692) = rxt(k,289)*y(k,45) + rxt(k,315)*y(k,49) + .500_r8*rxt(k,296) &
                      *y(k,51) + .650_r8*rxt(k,461)*y(k,184)
         mat(k,1163) = rxt(k,323)*y(k,130) + .300_r8*rxt(k,321)*y(k,204)  &
                      + .150_r8*rxt(k,322)*y(k,209)

         mat(k,1905) = -(rxt(k,182)*y(k,59) + (4._r8*rxt(k,259) + 4._r8*rxt(k,260) &
                      ) * y(k,204) + rxt(k,261)*y(k,209) + rxt(k,262)*y(k,130) &
                      + rxt(k,281)*y(k,200) + rxt(k,292)*y(k,203) + rxt(k,309) &
                      *y(k,201) + rxt(k,321)*y(k,226) + rxt(k,333)*y(k,215) + rxt(k,342) &
                      *y(k,217) + rxt(k,366)*y(k,211) + rxt(k,371)*y(k,212) + rxt(k,380) &
                      *y(k,107) + rxt(k,391)*y(k,231) + rxt(k,445)*y(k,221) + rxt(k,450) &
                      *y(k,227) + rxt(k,455)*y(k,228))
         mat(k,2062) = -rxt(k,182)*y(k,204)
         mat(k,2012) = -rxt(k,261)*y(k,204)
         mat(k,1854) = -rxt(k,262)*y(k,204)
         mat(k,830) = -rxt(k,281)*y(k,204)
         mat(k,1400) = -rxt(k,292)*y(k,204)
         mat(k,798) = -rxt(k,309)*y(k,204)
         mat(k,1166) = -rxt(k,321)*y(k,204)
         mat(k,1326) = -rxt(k,333)*y(k,204)
         mat(k,1368) = -rxt(k,342)*y(k,204)
         mat(k,1271) = -rxt(k,366)*y(k,204)
         mat(k,1304) = -rxt(k,371)*y(k,204)
         mat(k,1207) = -rxt(k,380)*y(k,204)
         mat(k,1184) = -rxt(k,391)*y(k,204)
         mat(k,1071) = -rxt(k,445)*y(k,204)
         mat(k,1130) = -rxt(k,450)*y(k,204)
         mat(k,1152) = -rxt(k,455)*y(k,204)

         mat(k,1019) = .280_r8*rxt(k,308)*y(k,140)
         mat(k,693) = rxt(k,295)*y(k,223)
         mat(k,450) = .700_r8*rxt(k,264)*y(k,223)
         mat(k,1444) = rxt(k,176)*y(k,56) + rxt(k,232)*y(k,79) + rxt(k,271)*y(k,222)  &
                      + rxt(k,265)*y(k,223)
         mat(k,2101) = rxt(k,176)*y(k,54)
         mat(k,886) = rxt(k,232)*y(k,54)
         mat(k,865) = .050_r8*rxt(k,385)*y(k,140)
         mat(k,1207) = mat(k,1207) + rxt(k,379)*y(k,203)
         mat(k,1854) = mat(k,1854) + rxt(k,294)*y(k,203) + .830_r8*rxt(k,411)*y(k,205)  &
                      + .170_r8*rxt(k,417)*y(k,218)
         mat(k,2163) = .280_r8*rxt(k,308)*y(k,29) + .050_r8*rxt(k,385)*y(k,104)
         mat(k,1400) = mat(k,1400) + rxt(k,379)*y(k,107) + rxt(k,294)*y(k,130)  &
                      + 4.000_r8*rxt(k,291)*y(k,203) + .900_r8*rxt(k,292)*y(k,204)  &
                      + .450_r8*rxt(k,293)*y(k,209) + rxt(k,365)*y(k,211) + rxt(k,370) &
                      *y(k,212) + rxt(k,332)*y(k,215) + rxt(k,341)*y(k,217)  &
                      + rxt(k,390)*y(k,231)
         mat(k,1905) = mat(k,1905) + .900_r8*rxt(k,292)*y(k,203)
         mat(k,770) = .830_r8*rxt(k,411)*y(k,130) + .330_r8*rxt(k,410)*y(k,209)
         mat(k,2012) = mat(k,2012) + .450_r8*rxt(k,293)*y(k,203) + .330_r8*rxt(k,410) &
                      *y(k,205) + .070_r8*rxt(k,416)*y(k,218)
         mat(k,1271) = mat(k,1271) + rxt(k,365)*y(k,203)
         mat(k,1304) = mat(k,1304) + rxt(k,370)*y(k,203)
         mat(k,1326) = mat(k,1326) + rxt(k,332)*y(k,203)
         mat(k,1368) = mat(k,1368) + rxt(k,341)*y(k,203)
         mat(k,875) = .170_r8*rxt(k,417)*y(k,130) + .070_r8*rxt(k,416)*y(k,209)
         mat(k,1539) = rxt(k,271)*y(k,54)
         mat(k,1703) = rxt(k,295)*y(k,50) + .700_r8*rxt(k,264)*y(k,53) + rxt(k,265) &
                      *y(k,54)
         mat(k,1184) = mat(k,1184) + rxt(k,390)*y(k,203)

         mat(k,765) = -(rxt(k,410)*y(k,209) + rxt(k,411)*y(k,130) + rxt(k,412) &
                      *y(k,131))
         mat(k,1968) = -rxt(k,410)*y(k,205)
         mat(k,1812) = -rxt(k,411)*y(k,205)
         mat(k,2208) = -rxt(k,412)*y(k,205)

         mat(k,570) = -((rxt(k,329) + rxt(k,330)) * y(k,130))
         mat(k,1801) = -(rxt(k,329) + rxt(k,330)) * y(k,206)

         mat(k,356) = rxt(k,328)*y(k,223)
         mat(k,1628) = rxt(k,328)*y(k,16)


         mat(k,1786) = .750_r8*rxt(k,298)*y(k,208)
         mat(k,719) = .750_r8*rxt(k,298)*y(k,130)

         mat(k,720) = -(rxt(k,297)*y(k,209) + rxt(k,298)*y(k,130))
         mat(k,1964) = -rxt(k,297)*y(k,208)
         mat(k,1808) = -rxt(k,298)*y(k,208)

         mat(k,555) = rxt(k,304)*y(k,223)
         mat(k,1645) = rxt(k,304)*y(k,25)

         mat(k,2013) = -((rxt(k,135) + rxt(k,136) + rxt(k,137)) * y(k,82) + rxt(k,139) &
                      *y(k,139) + rxt(k,140)*y(k,140) + rxt(k,144)*y(k,223) &
                      + 4._r8*rxt(k,149)*y(k,209) + rxt(k,159)*y(k,132) + rxt(k,164) &
                      *y(k,130) + rxt(k,169)*y(k,131) + (rxt(k,179) + rxt(k,180) &
                      ) * y(k,56) + rxt(k,186)*y(k,59) + rxt(k,212)*y(k,17) + rxt(k,218) &
                      *y(k,19) + rxt(k,255)*y(k,42) + rxt(k,261)*y(k,204) + rxt(k,268) &
                      *y(k,210) + rxt(k,282)*y(k,200) + rxt(k,293)*y(k,203) + rxt(k,297) &
                      *y(k,208) + rxt(k,310)*y(k,201) + rxt(k,318)*y(k,225) + rxt(k,322) &
                      *y(k,226) + rxt(k,334)*y(k,215) + rxt(k,343)*y(k,217) + rxt(k,347) &
                      *y(k,219) + rxt(k,357)*y(k,195) + rxt(k,367)*y(k,211) + rxt(k,372) &
                      *y(k,212) + rxt(k,381)*y(k,107) + rxt(k,392)*y(k,231) + rxt(k,396) &
                      *y(k,194) + rxt(k,399)*y(k,197) + rxt(k,403)*y(k,199) + rxt(k,406) &
                      *y(k,202) + rxt(k,410)*y(k,205) + rxt(k,413)*y(k,216) + rxt(k,416) &
                      *y(k,218) + rxt(k,419)*y(k,224) + rxt(k,426)*y(k,229) + rxt(k,432) &
                      *y(k,232) + rxt(k,435)*y(k,234) + rxt(k,446)*y(k,221) + rxt(k,451) &
                      *y(k,227) + rxt(k,456)*y(k,228))
         mat(k,1459) = -(rxt(k,135) + rxt(k,136) + rxt(k,137)) * y(k,209)
         mat(k,2262) = -rxt(k,139)*y(k,209)
         mat(k,2164) = -rxt(k,140)*y(k,209)
         mat(k,1704) = -rxt(k,144)*y(k,209)
         mat(k,1763) = -rxt(k,159)*y(k,209)
         mat(k,1855) = -rxt(k,164)*y(k,209)
         mat(k,2231) = -rxt(k,169)*y(k,209)
         mat(k,2102) = -(rxt(k,179) + rxt(k,180)) * y(k,209)
         mat(k,2063) = -rxt(k,186)*y(k,209)
         mat(k,1428) = -rxt(k,212)*y(k,209)
         mat(k,1499) = -rxt(k,218)*y(k,209)
         mat(k,2037) = -rxt(k,255)*y(k,209)
         mat(k,1906) = -rxt(k,261)*y(k,209)
         mat(k,444) = -rxt(k,268)*y(k,209)
         mat(k,831) = -rxt(k,282)*y(k,209)
         mat(k,1401) = -rxt(k,293)*y(k,209)
         mat(k,724) = -rxt(k,297)*y(k,209)
         mat(k,799) = -rxt(k,310)*y(k,209)
         mat(k,779) = -rxt(k,318)*y(k,209)
         mat(k,1167) = -rxt(k,322)*y(k,209)
         mat(k,1327) = -rxt(k,334)*y(k,209)
         mat(k,1369) = -rxt(k,343)*y(k,209)
         mat(k,688) = -rxt(k,347)*y(k,209)
         mat(k,928) = -rxt(k,357)*y(k,209)
         mat(k,1272) = -rxt(k,367)*y(k,209)
         mat(k,1305) = -rxt(k,372)*y(k,209)
         mat(k,1208) = -rxt(k,381)*y(k,209)
         mat(k,1185) = -rxt(k,392)*y(k,209)
         mat(k,528) = -rxt(k,396)*y(k,209)
         mat(k,502) = -rxt(k,399)*y(k,209)
         mat(k,439) = -rxt(k,403)*y(k,209)
         mat(k,634) = -rxt(k,406)*y(k,209)
         mat(k,771) = -rxt(k,410)*y(k,209)
         mat(k,731) = -rxt(k,413)*y(k,209)
         mat(k,876) = -rxt(k,416)*y(k,209)
         mat(k,464) = -rxt(k,419)*y(k,209)
         mat(k,746) = -rxt(k,426)*y(k,209)
         mat(k,763) = -rxt(k,432)*y(k,209)
         mat(k,510) = -rxt(k,435)*y(k,209)
         mat(k,1072) = -rxt(k,446)*y(k,209)
         mat(k,1131) = -rxt(k,451)*y(k,209)
         mat(k,1153) = -rxt(k,456)*y(k,209)

         mat(k,990) = .570_r8*rxt(k,440)*y(k,140)
         mat(k,173) = .650_r8*rxt(k,398)*y(k,223)
         mat(k,1428) = mat(k,1428) + rxt(k,211)*y(k,42)
         mat(k,1499) = mat(k,1499) + rxt(k,223)*y(k,223)
         mat(k,290) = .350_r8*rxt(k,277)*y(k,223)
         mat(k,559) = .130_r8*rxt(k,279)*y(k,140)
         mat(k,265) = rxt(k,284)*y(k,223)
         mat(k,1020) = .280_r8*rxt(k,308)*y(k,140)
         mat(k,2037) = mat(k,2037) + rxt(k,211)*y(k,17) + rxt(k,175)*y(k,56)  &
                      + rxt(k,256)*y(k,132) + rxt(k,257)*y(k,139)
         mat(k,596) = rxt(k,240)*y(k,56) + rxt(k,241)*y(k,223)
         mat(k,375) = rxt(k,243)*y(k,56) + rxt(k,244)*y(k,223)
         mat(k,110) = rxt(k,290)*y(k,223)
         mat(k,805) = rxt(k,263)*y(k,223)
         mat(k,1445) = rxt(k,272)*y(k,222)
         mat(k,2102) = mat(k,2102) + rxt(k,175)*y(k,42) + rxt(k,240)*y(k,43)  &
                      + rxt(k,243)*y(k,46) + rxt(k,178)*y(k,85)
         mat(k,2063) = mat(k,2063) + rxt(k,182)*y(k,204) + rxt(k,193)*y(k,223)
         mat(k,1090) = rxt(k,275)*y(k,223)
         mat(k,204) = .730_r8*rxt(k,409)*y(k,223)
         mat(k,308) = .500_r8*rxt(k,477)*y(k,223)
         mat(k,1056) = rxt(k,301)*y(k,223)
         mat(k,961) = rxt(k,302)*y(k,223)
         mat(k,602) = rxt(k,178)*y(k,56) + rxt(k,134)*y(k,139) + rxt(k,143)*y(k,223)
         mat(k,188) = rxt(k,266)*y(k,223)
         mat(k,947) = rxt(k,267)*y(k,223)
         mat(k,1105) = rxt(k,331)*y(k,223)
         mat(k,1114) = rxt(k,316)*y(k,223)
         mat(k,866) = .370_r8*rxt(k,385)*y(k,140)
         mat(k,584) = .300_r8*rxt(k,376)*y(k,223)
         mat(k,553) = rxt(k,377)*y(k,223)
         mat(k,1208) = mat(k,1208) + rxt(k,382)*y(k,130) + rxt(k,383)*y(k,132)  &
                      + rxt(k,379)*y(k,203) + 1.200_r8*rxt(k,380)*y(k,204)
         mat(k,395) = rxt(k,384)*y(k,223)
         mat(k,1220) = .140_r8*rxt(k,338)*y(k,140)
         mat(k,320) = .200_r8*rxt(k,340)*y(k,223)
         mat(k,613) = .500_r8*rxt(k,351)*y(k,223)
         mat(k,912) = .570_r8*rxt(k,443)*y(k,140)
         mat(k,1349) = .280_r8*rxt(k,352)*y(k,140)
         mat(k,384) = rxt(k,388)*y(k,223)
         mat(k,1038) = rxt(k,389)*y(k,223)
         mat(k,1855) = mat(k,1855) + rxt(k,382)*y(k,107) + rxt(k,358)*y(k,195)  &
                      + rxt(k,400)*y(k,197) + rxt(k,405)*y(k,199) + rxt(k,283) &
                      *y(k,200) + rxt(k,311)*y(k,201) + rxt(k,262)*y(k,204)  &
                      + .170_r8*rxt(k,411)*y(k,205) + rxt(k,329)*y(k,206)  &
                      + .250_r8*rxt(k,298)*y(k,208) + rxt(k,270)*y(k,210)  &
                      + .920_r8*rxt(k,368)*y(k,211) + .920_r8*rxt(k,374)*y(k,212)  &
                      + .470_r8*rxt(k,336)*y(k,215) + .400_r8*rxt(k,414)*y(k,216)  &
                      + .830_r8*rxt(k,417)*y(k,218) + rxt(k,420)*y(k,224) + rxt(k,319) &
                      *y(k,225) + .900_r8*rxt(k,452)*y(k,227) + .800_r8*rxt(k,457) &
                      *y(k,228) + rxt(k,427)*y(k,229) + rxt(k,393)*y(k,231)  &
                      + rxt(k,433)*y(k,232) + rxt(k,436)*y(k,234)
         mat(k,1763) = mat(k,1763) + rxt(k,256)*y(k,42) + rxt(k,383)*y(k,107)  &
                      + rxt(k,369)*y(k,211) + rxt(k,375)*y(k,212) + .470_r8*rxt(k,335) &
                      *y(k,215) + rxt(k,162)*y(k,223) + rxt(k,394)*y(k,231)
         mat(k,2262) = mat(k,2262) + rxt(k,257)*y(k,42) + rxt(k,134)*y(k,85)
         mat(k,2164) = mat(k,2164) + .570_r8*rxt(k,440)*y(k,6) + .130_r8*rxt(k,279) &
                      *y(k,25) + .280_r8*rxt(k,308)*y(k,29) + .370_r8*rxt(k,385) &
                      *y(k,104) + .140_r8*rxt(k,338)*y(k,111) + .570_r8*rxt(k,443) &
                      *y(k,116) + .280_r8*rxt(k,352)*y(k,117) + rxt(k,146)*y(k,223)
         mat(k,182) = .800_r8*rxt(k,421)*y(k,223)
         mat(k,847) = rxt(k,467)*y(k,223)
         mat(k,1083) = .200_r8*rxt(k,461)*y(k,223)
         mat(k,199) = .280_r8*rxt(k,429)*y(k,223)
         mat(k,221) = .380_r8*rxt(k,431)*y(k,223)
         mat(k,226) = .630_r8*rxt(k,437)*y(k,223)
         mat(k,928) = mat(k,928) + rxt(k,358)*y(k,130)
         mat(k,502) = mat(k,502) + rxt(k,400)*y(k,130)
         mat(k,439) = mat(k,439) + rxt(k,405)*y(k,130)
         mat(k,831) = mat(k,831) + rxt(k,283)*y(k,130) + 2.400_r8*rxt(k,280)*y(k,200)  &
                      + rxt(k,281)*y(k,204)
         mat(k,799) = mat(k,799) + rxt(k,311)*y(k,130) + rxt(k,309)*y(k,204)
         mat(k,1401) = mat(k,1401) + rxt(k,379)*y(k,107) + .900_r8*rxt(k,292)*y(k,204)  &
                      + rxt(k,365)*y(k,211) + rxt(k,370)*y(k,212) + .470_r8*rxt(k,332) &
                      *y(k,215) + rxt(k,390)*y(k,231)
         mat(k,1906) = mat(k,1906) + rxt(k,182)*y(k,59) + 1.200_r8*rxt(k,380)*y(k,107)  &
                      + rxt(k,262)*y(k,130) + rxt(k,281)*y(k,200) + rxt(k,309) &
                      *y(k,201) + .900_r8*rxt(k,292)*y(k,203) + 4.000_r8*rxt(k,259) &
                      *y(k,204) + rxt(k,366)*y(k,211) + rxt(k,371)*y(k,212)  &
                      + .730_r8*rxt(k,333)*y(k,215) + rxt(k,342)*y(k,217)  &
                      + .500_r8*rxt(k,445)*y(k,221) + .300_r8*rxt(k,321)*y(k,226)  &
                      + rxt(k,450)*y(k,227) + rxt(k,455)*y(k,228) + .800_r8*rxt(k,391) &
                      *y(k,231)
         mat(k,771) = mat(k,771) + .170_r8*rxt(k,411)*y(k,130) + .070_r8*rxt(k,410) &
                      *y(k,209)
         mat(k,575) = rxt(k,329)*y(k,130)
         mat(k,724) = mat(k,724) + .250_r8*rxt(k,298)*y(k,130)
         mat(k,2013) = mat(k,2013) + .070_r8*rxt(k,410)*y(k,205) + .160_r8*rxt(k,413) &
                      *y(k,216) + .330_r8*rxt(k,416)*y(k,218)
         mat(k,444) = mat(k,444) + rxt(k,270)*y(k,130)
         mat(k,1272) = mat(k,1272) + .920_r8*rxt(k,368)*y(k,130) + rxt(k,369)*y(k,132)  &
                      + rxt(k,365)*y(k,203) + rxt(k,366)*y(k,204)
         mat(k,1305) = mat(k,1305) + .920_r8*rxt(k,374)*y(k,130) + rxt(k,375)*y(k,132)  &
                      + rxt(k,370)*y(k,203) + rxt(k,371)*y(k,204)
         mat(k,1327) = mat(k,1327) + .470_r8*rxt(k,336)*y(k,130) + .470_r8*rxt(k,335) &
                      *y(k,132) + .470_r8*rxt(k,332)*y(k,203) + .730_r8*rxt(k,333) &
                      *y(k,204)
         mat(k,731) = mat(k,731) + .400_r8*rxt(k,414)*y(k,130) + .160_r8*rxt(k,413) &
                      *y(k,209)
         mat(k,1369) = mat(k,1369) + rxt(k,342)*y(k,204)
         mat(k,876) = mat(k,876) + .830_r8*rxt(k,417)*y(k,130) + .330_r8*rxt(k,416) &
                      *y(k,209)
         mat(k,1072) = mat(k,1072) + .500_r8*rxt(k,445)*y(k,204)
         mat(k,1540) = rxt(k,272)*y(k,54)
         mat(k,1704) = mat(k,1704) + .650_r8*rxt(k,398)*y(k,7) + rxt(k,223)*y(k,19)  &
                      + .350_r8*rxt(k,277)*y(k,24) + rxt(k,284)*y(k,26) + rxt(k,241) &
                      *y(k,43) + rxt(k,244)*y(k,46) + rxt(k,290)*y(k,47) + rxt(k,263) &
                      *y(k,52) + rxt(k,193)*y(k,59) + rxt(k,275)*y(k,62)  &
                      + .730_r8*rxt(k,409)*y(k,66) + .500_r8*rxt(k,477)*y(k,67)  &
                      + rxt(k,301)*y(k,80) + rxt(k,302)*y(k,81) + rxt(k,143)*y(k,85)  &
                      + rxt(k,266)*y(k,92) + rxt(k,267)*y(k,93) + rxt(k,331)*y(k,99)  &
                      + rxt(k,316)*y(k,101) + .300_r8*rxt(k,376)*y(k,105) + rxt(k,377) &
                      *y(k,106) + rxt(k,384)*y(k,108) + .200_r8*rxt(k,340)*y(k,112)  &
                      + .500_r8*rxt(k,351)*y(k,115) + rxt(k,388)*y(k,121) + rxt(k,389) &
                      *y(k,122) + rxt(k,162)*y(k,132) + rxt(k,146)*y(k,140)  &
                      + .800_r8*rxt(k,421)*y(k,148) + rxt(k,467)*y(k,157)  &
                      + .200_r8*rxt(k,461)*y(k,184) + .280_r8*rxt(k,429)*y(k,186)  &
                      + .380_r8*rxt(k,431)*y(k,188) + .630_r8*rxt(k,437)*y(k,190)
         mat(k,464) = mat(k,464) + rxt(k,420)*y(k,130)
         mat(k,779) = mat(k,779) + rxt(k,319)*y(k,130)
         mat(k,1167) = mat(k,1167) + .300_r8*rxt(k,321)*y(k,204)
         mat(k,1131) = mat(k,1131) + .900_r8*rxt(k,452)*y(k,130) + rxt(k,450)*y(k,204)
         mat(k,1153) = mat(k,1153) + .800_r8*rxt(k,457)*y(k,130) + rxt(k,455)*y(k,204)
         mat(k,746) = mat(k,746) + rxt(k,427)*y(k,130)
         mat(k,1185) = mat(k,1185) + rxt(k,393)*y(k,130) + rxt(k,394)*y(k,132)  &
                      + rxt(k,390)*y(k,203) + .800_r8*rxt(k,391)*y(k,204)
         mat(k,763) = mat(k,763) + rxt(k,433)*y(k,130)
         mat(k,510) = mat(k,510) + rxt(k,436)*y(k,130)

      end do

      end subroutine     nlnmat08

      subroutine     nlnmat09( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,441) = -(rxt(k,268)*y(k,209) + rxt(k,270)*y(k,130))
         mat(k,1944) = -rxt(k,268)*y(k,210)
         mat(k,1791) = -rxt(k,270)*y(k,210)

         mat(k,2022) = rxt(k,255)*y(k,209)
         mat(k,1944) = mat(k,1944) + rxt(k,255)*y(k,42)

         mat(k,1262) = -(rxt(k,365)*y(k,203) + rxt(k,366)*y(k,204) + rxt(k,367) &
                      *y(k,209) + rxt(k,368)*y(k,130) + rxt(k,369)*y(k,132))
         mat(k,1389) = -rxt(k,365)*y(k,211)
         mat(k,1892) = -rxt(k,366)*y(k,211)
         mat(k,1996) = -rxt(k,367)*y(k,211)
         mat(k,1839) = -rxt(k,368)*y(k,211)
         mat(k,1747) = -rxt(k,369)*y(k,211)

         mat(k,859) = .600_r8*rxt(k,386)*y(k,223)
         mat(k,1687) = .600_r8*rxt(k,386)*y(k,104)

         mat(k,1295) = -(rxt(k,370)*y(k,203) + rxt(k,371)*y(k,204) + rxt(k,372) &
                      *y(k,209) + rxt(k,374)*y(k,130) + rxt(k,375)*y(k,132))
         mat(k,1390) = -rxt(k,370)*y(k,212)
         mat(k,1893) = -rxt(k,371)*y(k,212)
         mat(k,1997) = -rxt(k,372)*y(k,212)
         mat(k,1840) = -rxt(k,374)*y(k,212)
         mat(k,1748) = -rxt(k,375)*y(k,212)

         mat(k,860) = .400_r8*rxt(k,386)*y(k,223)
         mat(k,1688) = .400_r8*rxt(k,386)*y(k,104)

         mat(k,72) = -(rxt(k,510)*y(k,209) + rxt(k,511)*y(k,130))
         mat(k,1922) = -rxt(k,510)*y(k,213)
         mat(k,1779) = -rxt(k,511)*y(k,213)

         mat(k,852) = rxt(k,513)*y(k,223)
         mat(k,1557) = rxt(k,513)*y(k,104)

         mat(k,78) = -(rxt(k,514)*y(k,209) + rxt(k,515)*y(k,130))
         mat(k,1923) = -rxt(k,514)*y(k,214)
         mat(k,1780) = -rxt(k,515)*y(k,214)

         mat(k,79) = rxt(k,516)*y(k,223)
         mat(k,1558) = rxt(k,516)*y(k,110)

         mat(k,1319) = -(rxt(k,332)*y(k,203) + rxt(k,333)*y(k,204) + rxt(k,334) &
                      *y(k,209) + rxt(k,335)*y(k,132) + (rxt(k,336) + rxt(k,337) &
                      ) * y(k,130))
         mat(k,1391) = -rxt(k,332)*y(k,215)
         mat(k,1894) = -rxt(k,333)*y(k,215)
         mat(k,1998) = -rxt(k,334)*y(k,215)
         mat(k,1749) = -rxt(k,335)*y(k,215)
         mat(k,1841) = -(rxt(k,336) + rxt(k,337)) * y(k,215)

         mat(k,1216) = .500_r8*rxt(k,339)*y(k,223)
         mat(k,317) = .200_r8*rxt(k,340)*y(k,223)
         mat(k,1339) = rxt(k,353)*y(k,223)
         mat(k,1689) = .500_r8*rxt(k,339)*y(k,111) + .200_r8*rxt(k,340)*y(k,112)  &
                      + rxt(k,353)*y(k,117)

         mat(k,727) = -(rxt(k,413)*y(k,209) + rxt(k,414)*y(k,130) + rxt(k,415) &
                      *y(k,131))
         mat(k,1965) = -rxt(k,413)*y(k,216)
         mat(k,1809) = -rxt(k,414)*y(k,216)
         mat(k,2207) = -rxt(k,415)*y(k,216)

         mat(k,1362) = -(rxt(k,341)*y(k,203) + rxt(k,342)*y(k,204) + rxt(k,343) &
                      *y(k,209) + 4._r8*rxt(k,344)*y(k,217) + rxt(k,345)*y(k,130) &
                      + rxt(k,346)*y(k,132) + rxt(k,354)*y(k,131))
         mat(k,1393) = -rxt(k,341)*y(k,217)
         mat(k,1896) = -rxt(k,342)*y(k,217)
         mat(k,2000) = -rxt(k,343)*y(k,217)
         mat(k,1843) = -rxt(k,345)*y(k,217)
         mat(k,1751) = -rxt(k,346)*y(k,217)
         mat(k,2219) = -rxt(k,354)*y(k,217)

         mat(k,1217) = .500_r8*rxt(k,339)*y(k,223)
         mat(k,318) = .500_r8*rxt(k,340)*y(k,223)
         mat(k,1691) = .500_r8*rxt(k,339)*y(k,111) + .500_r8*rxt(k,340)*y(k,112)

         mat(k,869) = -(rxt(k,416)*y(k,209) + rxt(k,417)*y(k,130) + rxt(k,418) &
                      *y(k,131))
         mat(k,1977) = -rxt(k,416)*y(k,218)
         mat(k,1818) = -rxt(k,417)*y(k,218)
         mat(k,2212) = -rxt(k,418)*y(k,218)

         mat(k,683) = -(rxt(k,347)*y(k,209) + rxt(k,348)*y(k,130))
         mat(k,1960) = -rxt(k,347)*y(k,219)
         mat(k,1807) = -rxt(k,348)*y(k,219)

         mat(k,513) = rxt(k,349)*y(k,223)
         mat(k,322) = rxt(k,350)*y(k,223)
         mat(k,1641) = rxt(k,349)*y(k,113) + rxt(k,350)*y(k,114)

         mat(k,86) = -(rxt(k,518)*y(k,209) + rxt(k,519)*y(k,130))
         mat(k,1924) = -rxt(k,518)*y(k,220)
         mat(k,1781) = -rxt(k,519)*y(k,220)

         mat(k,893) = rxt(k,521)*y(k,223)
         mat(k,1560) = rxt(k,521)*y(k,116)

         mat(k,1063) = -(rxt(k,445)*y(k,204) + rxt(k,446)*y(k,209) + rxt(k,447) &
                      *y(k,130) + rxt(k,448)*y(k,132))
         mat(k,1881) = -rxt(k,445)*y(k,221)
         mat(k,1985) = -rxt(k,446)*y(k,221)
         mat(k,1827) = -rxt(k,447)*y(k,221)
         mat(k,1734) = -rxt(k,448)*y(k,221)

         mat(k,979) = rxt(k,439)*y(k,132)
         mat(k,903) = rxt(k,442)*y(k,132)
         mat(k,1734) = mat(k,1734) + rxt(k,439)*y(k,6) + rxt(k,442)*y(k,116)  &
                      + .500_r8*rxt(k,459)*y(k,183)
         mat(k,387) = rxt(k,449)*y(k,223)
         mat(k,997) = .500_r8*rxt(k,459)*y(k,132)
         mat(k,1674) = rxt(k,449)*y(k,134)

         mat(k,1535) = -(rxt(k,125)*y(k,83) + rxt(k,126)*y(k,235) + rxt(k,129) &
                      *y(k,140) + (rxt(k,167) + rxt(k,168)) * y(k,119) + rxt(k,200) &
                      *y(k,33) + rxt(k,201)*y(k,34) + rxt(k,202)*y(k,36) + rxt(k,203) &
                      *y(k,37) + rxt(k,204)*y(k,38) + rxt(k,205)*y(k,39) + rxt(k,206) &
                      *y(k,40) + (rxt(k,207) + rxt(k,208)) * y(k,91) + rxt(k,227) &
                      *y(k,35) + rxt(k,228)*y(k,55) + rxt(k,229)*y(k,84) + (rxt(k,230) &
                      + rxt(k,231)) * y(k,87) + rxt(k,236)*y(k,64) + rxt(k,237) &
                      *y(k,65) + rxt(k,250)*y(k,41) + rxt(k,251)*y(k,43) + rxt(k,252) &
                      *y(k,88) + rxt(k,253)*y(k,89) + rxt(k,254)*y(k,90) + (rxt(k,271) &
                      + rxt(k,272) + rxt(k,273)) * y(k,54) + rxt(k,274)*y(k,92))
         mat(k,1413) = -rxt(k,125)*y(k,222)
         mat(k,2283) = -rxt(k,126)*y(k,222)
         mat(k,2159) = -rxt(k,129)*y(k,222)
         mat(k,190) = -(rxt(k,167) + rxt(k,168)) * y(k,222)
         mat(k,106) = -rxt(k,200)*y(k,222)
         mat(k,150) = -rxt(k,201)*y(k,222)
         mat(k,121) = -rxt(k,202)*y(k,222)
         mat(k,160) = -rxt(k,203)*y(k,222)
         mat(k,125) = -rxt(k,204)*y(k,222)
         mat(k,165) = -rxt(k,205)*y(k,222)
         mat(k,129) = -rxt(k,206)*y(k,222)
         mat(k,2182) = -(rxt(k,207) + rxt(k,208)) * y(k,222)
         mat(k,156) = -rxt(k,227)*y(k,222)
         mat(k,455) = -rxt(k,228)*y(k,222)
         mat(k,114) = -rxt(k,229)*y(k,222)
         mat(k,819) = -(rxt(k,230) + rxt(k,231)) * y(k,222)
         mat(k,255) = -rxt(k,236)*y(k,222)
         mat(k,234) = -rxt(k,237)*y(k,222)
         mat(k,475) = -rxt(k,250)*y(k,222)
         mat(k,593) = -rxt(k,251)*y(k,222)
         mat(k,229) = -rxt(k,252)*y(k,222)
         mat(k,259) = -rxt(k,253)*y(k,222)
         mat(k,312) = -rxt(k,254)*y(k,222)
         mat(k,1441) = -(rxt(k,271) + rxt(k,272) + rxt(k,273)) * y(k,222)
         mat(k,186) = -rxt(k,274)*y(k,222)

         mat(k,1700) = -(rxt(k,142)*y(k,83) + rxt(k,143)*y(k,85) + rxt(k,144)*y(k,209) &
                      + rxt(k,145)*y(k,139) + rxt(k,146)*y(k,140) + (4._r8*rxt(k,147) &
                      + 4._r8*rxt(k,148)) * y(k,223) + rxt(k,150)*y(k,96) + rxt(k,162) &
                      *y(k,132) + rxt(k,163)*y(k,118) + rxt(k,171)*y(k,131) + rxt(k,172) &
                      *y(k,95) + rxt(k,191)*y(k,60) + (rxt(k,193) + rxt(k,194) &
                      ) * y(k,59) + rxt(k,196)*y(k,91) + rxt(k,199)*y(k,98) + rxt(k,223) &
                      *y(k,19) + rxt(k,225)*y(k,87) + rxt(k,239)*y(k,41) + rxt(k,241) &
                      *y(k,43) + rxt(k,242)*y(k,44) + rxt(k,244)*y(k,46) + rxt(k,246) &
                      *y(k,55) + rxt(k,247)*y(k,88) + rxt(k,248)*y(k,89) + rxt(k,249) &
                      *y(k,90) + rxt(k,258)*y(k,42) + rxt(k,263)*y(k,52) + rxt(k,264) &
                      *y(k,53) + rxt(k,265)*y(k,54) + rxt(k,266)*y(k,92) + rxt(k,267) &
                      *y(k,93) + rxt(k,275)*y(k,62) + rxt(k,277)*y(k,24) + rxt(k,284) &
                      *y(k,26) + rxt(k,285)*y(k,27) + rxt(k,287)*y(k,28) + rxt(k,289) &
                      *y(k,45) + rxt(k,290)*y(k,47) + rxt(k,295)*y(k,50) + rxt(k,296) &
                      *y(k,51) + rxt(k,301)*y(k,80) + rxt(k,302)*y(k,81) + rxt(k,303) &
                      *y(k,145) + rxt(k,304)*y(k,25) + rxt(k,312)*y(k,30) + rxt(k,313) &
                      *y(k,31) + rxt(k,315)*y(k,49) + rxt(k,316)*y(k,101) + rxt(k,317) &
                      *y(k,133) + rxt(k,320)*y(k,152) + rxt(k,324)*y(k,153) + rxt(k,325) &
                      *y(k,29) + rxt(k,326)*y(k,48) + rxt(k,328)*y(k,16) + rxt(k,331) &
                      *y(k,99) + rxt(k,339)*y(k,111) + rxt(k,340)*y(k,112) + rxt(k,349) &
                      *y(k,113) + rxt(k,350)*y(k,114) + rxt(k,351)*y(k,115) + rxt(k,353) &
                      *y(k,117) + rxt(k,356)*y(k,1) + rxt(k,360)*y(k,2) + rxt(k,361) &
                      *y(k,15) + rxt(k,362)*y(k,100) + rxt(k,363)*y(k,102) + rxt(k,364) &
                      *y(k,103) + rxt(k,376)*y(k,105) + rxt(k,377)*y(k,106) + rxt(k,384) &
                      *y(k,108) + rxt(k,386)*y(k,104) + rxt(k,387)*y(k,109) + rxt(k,388) &
                      *y(k,121) + rxt(k,389)*y(k,122) + rxt(k,395)*y(k,187) + rxt(k,398) &
                      *y(k,7) + rxt(k,401)*y(k,8) + rxt(k,402)*y(k,22) + rxt(k,404) &
                      *y(k,23) + rxt(k,408)*y(k,32) + rxt(k,409)*y(k,66) + rxt(k,421) &
                      *y(k,148) + rxt(k,424)*y(k,149) + rxt(k,428)*y(k,185) + rxt(k,429) &
                      *y(k,186) + rxt(k,431)*y(k,188) + rxt(k,434)*y(k,189) + rxt(k,437) &
                      *y(k,190) + rxt(k,438)*y(k,191) + rxt(k,441)*y(k,6) + rxt(k,444) &
                      *y(k,116) + rxt(k,449)*y(k,134) + rxt(k,453)*y(k,180) + rxt(k,454) &
                      *y(k,181) + rxt(k,458)*y(k,182) + rxt(k,460)*y(k,183) + rxt(k,461) &
                      *y(k,184) + (rxt(k,463) + rxt(k,477)) * y(k,67) + rxt(k,465) &
                      *y(k,143) + rxt(k,467)*y(k,157) + rxt(k,471)*y(k,154) + rxt(k,476) &
                      *y(k,156) + rxt(k,479)*y(k,126))
         mat(k,1414) = -rxt(k,142)*y(k,223)
         mat(k,601) = -rxt(k,143)*y(k,223)
         mat(k,2009) = -rxt(k,144)*y(k,223)
         mat(k,2258) = -rxt(k,145)*y(k,223)
         mat(k,2160) = -rxt(k,146)*y(k,223)
         mat(k,422) = -rxt(k,150)*y(k,223)
         mat(k,1759) = -rxt(k,162)*y(k,223)
         mat(k,484) = -rxt(k,163)*y(k,223)
         mat(k,2227) = -rxt(k,171)*y(k,223)
         mat(k,1473) = -rxt(k,172)*y(k,223)
         mat(k,935) = -rxt(k,191)*y(k,223)
         mat(k,2059) = -(rxt(k,193) + rxt(k,194)) * y(k,223)
         mat(k,2183) = -rxt(k,196)*y(k,223)
         mat(k,809) = -rxt(k,199)*y(k,223)
         mat(k,1495) = -rxt(k,223)*y(k,223)
         mat(k,820) = -rxt(k,225)*y(k,223)
         mat(k,476) = -rxt(k,239)*y(k,223)
         mat(k,594) = -rxt(k,241)*y(k,223)
         mat(k,132) = -rxt(k,242)*y(k,223)
         mat(k,373) = -rxt(k,244)*y(k,223)
         mat(k,456) = -rxt(k,246)*y(k,223)
         mat(k,230) = -rxt(k,247)*y(k,223)
         mat(k,260) = -rxt(k,248)*y(k,223)
         mat(k,313) = -rxt(k,249)*y(k,223)
         mat(k,2033) = -rxt(k,258)*y(k,223)
         mat(k,804) = -rxt(k,263)*y(k,223)
         mat(k,449) = -rxt(k,264)*y(k,223)
         mat(k,1442) = -rxt(k,265)*y(k,223)
         mat(k,187) = -rxt(k,266)*y(k,223)
         mat(k,946) = -rxt(k,267)*y(k,223)
         mat(k,1089) = -rxt(k,275)*y(k,223)
         mat(k,289) = -rxt(k,277)*y(k,223)
         mat(k,264) = -rxt(k,284)*y(k,223)
         mat(k,348) = -rxt(k,285)*y(k,223)
         mat(k,296) = -rxt(k,287)*y(k,223)
         mat(k,1046) = -rxt(k,289)*y(k,223)
         mat(k,109) = -rxt(k,290)*y(k,223)
         mat(k,692) = -rxt(k,295)*y(k,223)
         mat(k,619) = -rxt(k,296)*y(k,223)
         mat(k,1055) = -rxt(k,301)*y(k,223)
         mat(k,960) = -rxt(k,302)*y(k,223)
         mat(k,533) = -rxt(k,303)*y(k,223)
         mat(k,558) = -rxt(k,304)*y(k,223)
         mat(k,412) = -rxt(k,312)*y(k,223)
         mat(k,117) = -rxt(k,313)*y(k,223)
         mat(k,1229) = -rxt(k,315)*y(k,223)
         mat(k,1113) = -rxt(k,316)*y(k,223)
         mat(k,840) = -rxt(k,317)*y(k,223)
         mat(k,542) = -rxt(k,320)*y(k,223)
         mat(k,406) = -rxt(k,324)*y(k,223)
         mat(k,1016) = -rxt(k,325)*y(k,223)
         mat(k,953) = -rxt(k,326)*y(k,223)
         mat(k,359) = -rxt(k,328)*y(k,223)
         mat(k,1102) = -rxt(k,331)*y(k,223)
         mat(k,1219) = -rxt(k,339)*y(k,223)
         mat(k,319) = -rxt(k,340)*y(k,223)
         mat(k,516) = -rxt(k,349)*y(k,223)
         mat(k,325) = -rxt(k,350)*y(k,223)
         mat(k,611) = -rxt(k,351)*y(k,223)
         mat(k,1345) = -rxt(k,353)*y(k,223)
         mat(k,644) = -rxt(k,356)*y(k,223)
         mat(k,679) = -rxt(k,360)*y(k,223)
         mat(k,252) = -rxt(k,361)*y(k,223)
         mat(k,239) = -rxt(k,362)*y(k,223)
         mat(k,334) = -rxt(k,363)*y(k,223)
         mat(k,143) = -rxt(k,364)*y(k,223)
         mat(k,583) = -rxt(k,376)*y(k,223)
         mat(k,552) = -rxt(k,377)*y(k,223)
         mat(k,394) = -rxt(k,384)*y(k,223)
         mat(k,863) = -rxt(k,386)*y(k,223)
         mat(k,700) = -rxt(k,387)*y(k,223)
         mat(k,383) = -rxt(k,388)*y(k,223)
         mat(k,1035) = -rxt(k,389)*y(k,223)
         mat(k,211) = -rxt(k,395)*y(k,223)
         mat(k,172) = -rxt(k,398)*y(k,223)
         mat(k,401) = -rxt(k,401)*y(k,223)
         mat(k,246) = -rxt(k,402)*y(k,223)
         mat(k,330) = -rxt(k,404)*y(k,223)
         mat(k,274) = -rxt(k,408)*y(k,223)
         mat(k,203) = -rxt(k,409)*y(k,223)
         mat(k,181) = -rxt(k,421)*y(k,223)
         mat(k,343) = -rxt(k,424)*y(k,223)
         mat(k,669) = -rxt(k,428)*y(k,223)
         mat(k,198) = -rxt(k,429)*y(k,223)
         mat(k,220) = -rxt(k,431)*y(k,223)
         mat(k,716) = -rxt(k,434)*y(k,223)
         mat(k,225) = -rxt(k,437)*y(k,223)
         mat(k,431) = -rxt(k,438)*y(k,223)
         mat(k,986) = -rxt(k,441)*y(k,223)
         mat(k,910) = -rxt(k,444)*y(k,223)
         mat(k,389) = -rxt(k,449)*y(k,223)
         mat(k,655) = -rxt(k,453)*y(k,223)
         mat(k,626) = -rxt(k,454)*y(k,223)
         mat(k,492) = -rxt(k,458)*y(k,223)
         mat(k,1001) = -rxt(k,460)*y(k,223)
         mat(k,1081) = -rxt(k,461)*y(k,223)
         mat(k,306) = -(rxt(k,463) + rxt(k,477)) * y(k,223)
         mat(k,369) = -rxt(k,465)*y(k,223)
         mat(k,846) = -rxt(k,467)*y(k,223)
         mat(k,520) = -rxt(k,471)*y(k,223)
         mat(k,1241) = -rxt(k,476)*y(k,223)
         mat(k,103) = -rxt(k,479)*y(k,223)

         mat(k,986) = mat(k,986) + .630_r8*rxt(k,440)*y(k,140)
         mat(k,289) = mat(k,289) + .650_r8*rxt(k,277)*y(k,223)
         mat(k,558) = mat(k,558) + .130_r8*rxt(k,279)*y(k,140)
         mat(k,348) = mat(k,348) + .500_r8*rxt(k,285)*y(k,223)
         mat(k,1016) = mat(k,1016) + .360_r8*rxt(k,308)*y(k,140)
         mat(k,2033) = mat(k,2033) + rxt(k,257)*y(k,139)
         mat(k,449) = mat(k,449) + .300_r8*rxt(k,264)*y(k,223)
         mat(k,1442) = mat(k,1442) + rxt(k,271)*y(k,222)
         mat(k,2098) = rxt(k,180)*y(k,209)
         mat(k,884) = rxt(k,234)*y(k,235)
         mat(k,1456) = rxt(k,141)*y(k,140) + 2.000_r8*rxt(k,136)*y(k,209)
         mat(k,1414) = mat(k,1414) + rxt(k,133)*y(k,139) + rxt(k,125)*y(k,222)
         mat(k,601) = mat(k,601) + rxt(k,134)*y(k,139)
         mat(k,820) = mat(k,820) + rxt(k,224)*y(k,139) + rxt(k,230)*y(k,222)
         mat(k,2183) = mat(k,2183) + rxt(k,195)*y(k,139) + rxt(k,207)*y(k,222)
         mat(k,187) = mat(k,187) + rxt(k,274)*y(k,222)
         mat(k,787) = rxt(k,226)*y(k,139)
         mat(k,809) = mat(k,809) + rxt(k,198)*y(k,139)
         mat(k,863) = mat(k,863) + .320_r8*rxt(k,385)*y(k,140)
         mat(k,700) = mat(k,700) + .600_r8*rxt(k,387)*y(k,223)
         mat(k,1219) = mat(k,1219) + .240_r8*rxt(k,338)*y(k,140)
         mat(k,319) = mat(k,319) + .100_r8*rxt(k,340)*y(k,223)
         mat(k,910) = mat(k,910) + .630_r8*rxt(k,443)*y(k,140)
         mat(k,1345) = mat(k,1345) + .360_r8*rxt(k,352)*y(k,140)
         mat(k,1851) = rxt(k,164)*y(k,209)
         mat(k,1759) = mat(k,1759) + rxt(k,159)*y(k,209)
         mat(k,2258) = mat(k,2258) + rxt(k,257)*y(k,42) + rxt(k,133)*y(k,83)  &
                      + rxt(k,134)*y(k,85) + rxt(k,224)*y(k,87) + rxt(k,195)*y(k,91)  &
                      + rxt(k,226)*y(k,97) + rxt(k,198)*y(k,98) + rxt(k,139)*y(k,209)
         mat(k,2160) = mat(k,2160) + .630_r8*rxt(k,440)*y(k,6) + .130_r8*rxt(k,279) &
                      *y(k,25) + .360_r8*rxt(k,308)*y(k,29) + rxt(k,141)*y(k,82)  &
                      + .320_r8*rxt(k,385)*y(k,104) + .240_r8*rxt(k,338)*y(k,111)  &
                      + .630_r8*rxt(k,443)*y(k,116) + .360_r8*rxt(k,352)*y(k,117)  &
                      + rxt(k,140)*y(k,209)
         mat(k,542) = mat(k,542) + .500_r8*rxt(k,320)*y(k,223)
         mat(k,211) = mat(k,211) + .500_r8*rxt(k,395)*y(k,223)
         mat(k,526) = .400_r8*rxt(k,396)*y(k,209)
         mat(k,1397) = .450_r8*rxt(k,293)*y(k,209)
         mat(k,768) = .400_r8*rxt(k,410)*y(k,209)
         mat(k,2009) = mat(k,2009) + rxt(k,180)*y(k,56) + 2.000_r8*rxt(k,136)*y(k,82)  &
                      + rxt(k,164)*y(k,130) + rxt(k,159)*y(k,132) + rxt(k,139) &
                      *y(k,139) + rxt(k,140)*y(k,140) + .400_r8*rxt(k,396)*y(k,194)  &
                      + .450_r8*rxt(k,293)*y(k,203) + .400_r8*rxt(k,410)*y(k,205)  &
                      + .450_r8*rxt(k,343)*y(k,217) + .400_r8*rxt(k,416)*y(k,218)  &
                      + .200_r8*rxt(k,347)*y(k,219) + .150_r8*rxt(k,322)*y(k,226)
         mat(k,1365) = .450_r8*rxt(k,343)*y(k,209)
         mat(k,873) = .400_r8*rxt(k,416)*y(k,209)
         mat(k,686) = .200_r8*rxt(k,347)*y(k,209)
         mat(k,1536) = rxt(k,271)*y(k,54) + rxt(k,125)*y(k,83) + rxt(k,230)*y(k,87)  &
                      + rxt(k,207)*y(k,91) + rxt(k,274)*y(k,92) + 2.000_r8*rxt(k,126) &
                      *y(k,235)
         mat(k,1700) = mat(k,1700) + .650_r8*rxt(k,277)*y(k,24) + .500_r8*rxt(k,285) &
                      *y(k,27) + .300_r8*rxt(k,264)*y(k,53) + .600_r8*rxt(k,387) &
                      *y(k,109) + .100_r8*rxt(k,340)*y(k,112) + .500_r8*rxt(k,320) &
                      *y(k,152) + .500_r8*rxt(k,395)*y(k,187)
         mat(k,1164) = .150_r8*rxt(k,322)*y(k,209)
         mat(k,2284) = rxt(k,234)*y(k,79) + 2.000_r8*rxt(k,126)*y(k,222)

      end do

      end subroutine     nlnmat09

      subroutine     nlnmat10( avec_len, mat, y, rxt )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,460) = -(rxt(k,419)*y(k,209) + rxt(k,420)*y(k,130))
         mat(k,1946) = -rxt(k,419)*y(k,224)
         mat(k,1792) = -rxt(k,420)*y(k,224)

         mat(k,201) = .200_r8*rxt(k,409)*y(k,223)
         mat(k,179) = .140_r8*rxt(k,421)*y(k,223)
         mat(k,341) = rxt(k,424)*y(k,223)
         mat(k,1613) = .200_r8*rxt(k,409)*y(k,66) + .140_r8*rxt(k,421)*y(k,148)  &
                      + rxt(k,424)*y(k,149)

         mat(k,774) = -(rxt(k,318)*y(k,209) + rxt(k,319)*y(k,130))
         mat(k,1969) = -rxt(k,318)*y(k,225)
         mat(k,1813) = -rxt(k,319)*y(k,225)

         mat(k,1005) = rxt(k,325)*y(k,223)
         mat(k,539) = .500_r8*rxt(k,320)*y(k,223)
         mat(k,1650) = rxt(k,325)*y(k,29) + .500_r8*rxt(k,320)*y(k,152)

         mat(k,1161) = -(rxt(k,321)*y(k,204) + rxt(k,322)*y(k,209) + rxt(k,323) &
                      *y(k,130))
         mat(k,1887) = -rxt(k,321)*y(k,226)
         mat(k,1991) = -rxt(k,322)*y(k,226)
         mat(k,1834) = -rxt(k,323)*y(k,226)

         mat(k,983) = .060_r8*rxt(k,440)*y(k,140)
         mat(k,951) = rxt(k,326)*y(k,223)
         mat(k,907) = .060_r8*rxt(k,443)*y(k,140)
         mat(k,2143) = .060_r8*rxt(k,440)*y(k,6) + .060_r8*rxt(k,443)*y(k,116)
         mat(k,404) = rxt(k,324)*y(k,223)
         mat(k,1079) = .150_r8*rxt(k,461)*y(k,223)
         mat(k,1681) = rxt(k,326)*y(k,48) + rxt(k,324)*y(k,153) + .150_r8*rxt(k,461) &
                      *y(k,184)

         mat(k,1124) = -(rxt(k,450)*y(k,204) + rxt(k,451)*y(k,209) + rxt(k,452) &
                      *y(k,130))
         mat(k,1885) = -rxt(k,450)*y(k,227)
         mat(k,1989) = -rxt(k,451)*y(k,227)
         mat(k,1832) = -rxt(k,452)*y(k,227)

         mat(k,1739) = .500_r8*rxt(k,459)*y(k,183)
         mat(k,654) = rxt(k,453)*y(k,223)
         mat(k,1000) = .500_r8*rxt(k,459)*y(k,132) + rxt(k,460)*y(k,223)
         mat(k,1679) = rxt(k,453)*y(k,180) + rxt(k,460)*y(k,183)

         mat(k,1145) = -(rxt(k,455)*y(k,204) + rxt(k,456)*y(k,209) + rxt(k,457) &
                      *y(k,130))
         mat(k,1886) = -rxt(k,455)*y(k,228)
         mat(k,1990) = -rxt(k,456)*y(k,228)
         mat(k,1833) = -rxt(k,457)*y(k,228)

         mat(k,982) = rxt(k,441)*y(k,223)
         mat(k,906) = rxt(k,444)*y(k,223)
         mat(k,491) = rxt(k,458)*y(k,223)
         mat(k,1680) = rxt(k,441)*y(k,6) + rxt(k,444)*y(k,116) + rxt(k,458)*y(k,182)

         mat(k,738) = -(rxt(k,426)*y(k,209) + rxt(k,427)*y(k,130))
         mat(k,1966) = -rxt(k,426)*y(k,229)
         mat(k,1810) = -rxt(k,427)*y(k,229)

         mat(k,663) = rxt(k,428)*y(k,223)
         mat(k,197) = .650_r8*rxt(k,429)*y(k,223)
         mat(k,1647) = rxt(k,428)*y(k,185) + .650_r8*rxt(k,429)*y(k,186)

         mat(k,92) = -(rxt(k,524)*y(k,209) + rxt(k,525)*y(k,130))
         mat(k,1925) = -rxt(k,524)*y(k,230)
         mat(k,1782) = -rxt(k,525)*y(k,230)

         mat(k,192) = rxt(k,523)*y(k,223)
         mat(k,1561) = rxt(k,523)*y(k,186)

         mat(k,1178) = -(rxt(k,390)*y(k,203) + rxt(k,391)*y(k,204) + rxt(k,392) &
                      *y(k,209) + rxt(k,393)*y(k,130) + rxt(k,394)*y(k,132))
         mat(k,1385) = -rxt(k,390)*y(k,231)
         mat(k,1888) = -rxt(k,391)*y(k,231)
         mat(k,1992) = -rxt(k,392)*y(k,231)
         mat(k,1835) = -rxt(k,393)*y(k,231)
         mat(k,1742) = -rxt(k,394)*y(k,231)

         mat(k,238) = rxt(k,362)*y(k,223)
         mat(k,333) = rxt(k,363)*y(k,223)
         mat(k,142) = rxt(k,364)*y(k,223)
         mat(k,697) = .400_r8*rxt(k,387)*y(k,223)
         mat(k,210) = .500_r8*rxt(k,395)*y(k,223)
         mat(k,1682) = rxt(k,362)*y(k,100) + rxt(k,363)*y(k,102) + rxt(k,364)*y(k,103)  &
                      + .400_r8*rxt(k,387)*y(k,109) + .500_r8*rxt(k,395)*y(k,187)

         mat(k,754) = -(rxt(k,432)*y(k,209) + rxt(k,433)*y(k,130))
         mat(k,1967) = -rxt(k,432)*y(k,232)
         mat(k,1811) = -rxt(k,433)*y(k,232)

         mat(k,217) = .560_r8*rxt(k,431)*y(k,223)
         mat(k,709) = rxt(k,434)*y(k,223)
         mat(k,1648) = .560_r8*rxt(k,431)*y(k,188) + rxt(k,434)*y(k,189)

         mat(k,98) = -(rxt(k,527)*y(k,209) + rxt(k,528)*y(k,130))
         mat(k,1926) = -rxt(k,527)*y(k,233)
         mat(k,1783) = -rxt(k,528)*y(k,233)

         mat(k,212) = rxt(k,526)*y(k,223)
         mat(k,1562) = rxt(k,526)*y(k,188)

         mat(k,505) = -(rxt(k,435)*y(k,209) + rxt(k,436)*y(k,130))
         mat(k,1951) = -rxt(k,435)*y(k,234)
         mat(k,1797) = -rxt(k,436)*y(k,234)

         mat(k,224) = .300_r8*rxt(k,437)*y(k,223)
         mat(k,428) = rxt(k,438)*y(k,223)
         mat(k,1620) = .300_r8*rxt(k,437)*y(k,190) + rxt(k,438)*y(k,191)

         mat(k,2296) = -(rxt(k,126)*y(k,222) + rxt(k,234)*y(k,79) + rxt(k,478) &
                      *y(k,158))
         mat(k,1548) = -rxt(k,126)*y(k,235)
         mat(k,887) = -rxt(k,234)*y(k,235)
         mat(k,270) = -rxt(k,478)*y(k,235)

         mat(k,299) = rxt(k,287)*y(k,223)
         mat(k,414) = rxt(k,312)*y(k,223)
         mat(k,118) = rxt(k,313)*y(k,223)
         mat(k,479) = rxt(k,239)*y(k,223)
         mat(k,2045) = rxt(k,258)*y(k,223)
         mat(k,599) = rxt(k,241)*y(k,223)
         mat(k,134) = rxt(k,242)*y(k,223)
         mat(k,1050) = rxt(k,289)*y(k,223)
         mat(k,378) = rxt(k,244)*y(k,223)
         mat(k,955) = rxt(k,326)*y(k,223)
         mat(k,1233) = rxt(k,315)*y(k,223)
         mat(k,694) = rxt(k,295)*y(k,223)
         mat(k,622) = rxt(k,296)*y(k,223)
         mat(k,452) = rxt(k,264)*y(k,223)
         mat(k,1450) = rxt(k,265)*y(k,223)
         mat(k,1465) = rxt(k,137)*y(k,209)
         mat(k,1420) = rxt(k,142)*y(k,223)
         mat(k,606) = rxt(k,143)*y(k,223)
         mat(k,822) = rxt(k,225)*y(k,223)
         mat(k,315) = rxt(k,249)*y(k,223)
         mat(k,2195) = (rxt(k,537)+rxt(k,542))*y(k,97) + (rxt(k,530)+rxt(k,536) &
                       +rxt(k,541))*y(k,98) + rxt(k,196)*y(k,223)
         mat(k,949) = rxt(k,267)*y(k,223)
         mat(k,1483) = rxt(k,172)*y(k,223)
         mat(k,426) = rxt(k,150)*y(k,223)
         mat(k,791) = (rxt(k,537)+rxt(k,542))*y(k,91)
         mat(k,814) = (rxt(k,530)+rxt(k,536)+rxt(k,541))*y(k,91) + rxt(k,199)*y(k,223)
         mat(k,1224) = .500_r8*rxt(k,339)*y(k,223)
         mat(k,104) = rxt(k,479)*y(k,223)
         mat(k,545) = rxt(k,320)*y(k,223)
         mat(k,408) = rxt(k,324)*y(k,223)
         mat(k,2021) = rxt(k,137)*y(k,82) + rxt(k,144)*y(k,223)
         mat(k,1712) = rxt(k,287)*y(k,28) + rxt(k,312)*y(k,30) + rxt(k,313)*y(k,31)  &
                      + rxt(k,239)*y(k,41) + rxt(k,258)*y(k,42) + rxt(k,241)*y(k,43)  &
                      + rxt(k,242)*y(k,44) + rxt(k,289)*y(k,45) + rxt(k,244)*y(k,46)  &
                      + rxt(k,326)*y(k,48) + rxt(k,315)*y(k,49) + rxt(k,295)*y(k,50)  &
                      + rxt(k,296)*y(k,51) + rxt(k,264)*y(k,53) + rxt(k,265)*y(k,54)  &
                      + rxt(k,142)*y(k,83) + rxt(k,143)*y(k,85) + rxt(k,225)*y(k,87)  &
                      + rxt(k,249)*y(k,90) + rxt(k,196)*y(k,91) + rxt(k,267)*y(k,93)  &
                      + rxt(k,172)*y(k,95) + rxt(k,150)*y(k,96) + rxt(k,199)*y(k,98)  &
                      + .500_r8*rxt(k,339)*y(k,111) + rxt(k,479)*y(k,126) + rxt(k,320) &
                      *y(k,152) + rxt(k,324)*y(k,153) + rxt(k,144)*y(k,209)  &
                      + 2.000_r8*rxt(k,147)*y(k,223)

      end do

      end subroutine     nlnmat10

      subroutine     nlnmat_finit( avec_len, mat, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  dti(veclen)
      real(r8), intent(in)    ::  lmat(veclen,nzcnt)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)


!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k

!----------------------------------------------
!       ... complete matrix entries implicit species
!----------------------------------------------

      do k = 1,avec_len
         mat(k,   1) = lmat(k,   1)
         mat(k,   2) = lmat(k,   2)
         mat(k,   3) = lmat(k,   3)
         mat(k,   4) = lmat(k,   4)
         mat(k,   5) = lmat(k,   5)
         mat(k,   6) = lmat(k,   6)
         mat(k,   7) = lmat(k,   7)
         mat(k,   8) = lmat(k,   8)
         mat(k,   9) = lmat(k,   9)
         mat(k,  10) = lmat(k,  10)
         mat(k,  11) = lmat(k,  11)
         mat(k,  12) = lmat(k,  12)
         mat(k,  13) = lmat(k,  13)
         mat(k,  14) = lmat(k,  14)
         mat(k,  15) = lmat(k,  15)
         mat(k,  16) = lmat(k,  16)
         mat(k,  17) = lmat(k,  17)
         mat(k,  18) = lmat(k,  18)
         mat(k,  19) = lmat(k,  19)
         mat(k,  20) = lmat(k,  20)
         mat(k,  21) = lmat(k,  21)
         mat(k,  22) = lmat(k,  22)
         mat(k,  23) = lmat(k,  23)
         mat(k,  24) = lmat(k,  24)
         mat(k,  25) = lmat(k,  25)
         mat(k,  26) = lmat(k,  26)
         mat(k,  27) = lmat(k,  27)
         mat(k,  28) = lmat(k,  28)
         mat(k,  29) = lmat(k,  29)
         mat(k,  30) = lmat(k,  30)
         mat(k,  31) = lmat(k,  31)
         mat(k,  32) = lmat(k,  32)
         mat(k,  33) = lmat(k,  33)
         mat(k,  34) = lmat(k,  34)
         mat(k,  35) = lmat(k,  35)
         mat(k,  36) = lmat(k,  36)
         mat(k,  37) = lmat(k,  37)
         mat(k,  38) = lmat(k,  38)
         mat(k,  39) = lmat(k,  39)
         mat(k,  40) = lmat(k,  40)
         mat(k,  41) = lmat(k,  41)
         mat(k,  42) = lmat(k,  42)
         mat(k,  43) = lmat(k,  43)
         mat(k,  44) = lmat(k,  44)
         mat(k,  45) = lmat(k,  45)
         mat(k,  46) = lmat(k,  46)
         mat(k,  47) = lmat(k,  47)
         mat(k,  48) = lmat(k,  48)
         mat(k,  54) = mat(k,  54) + lmat(k,  54)
         mat(k,  60) = mat(k,  60) + lmat(k,  60)
         mat(k,  66) = mat(k,  66) + lmat(k,  66)
         mat(k,  72) = mat(k,  72) + lmat(k,  72)
         mat(k,  78) = mat(k,  78) + lmat(k,  78)
         mat(k,  80) = mat(k,  80) + lmat(k,  80)
         mat(k,  86) = mat(k,  86) + lmat(k,  86)
         mat(k,  92) = mat(k,  92) + lmat(k,  92)
         mat(k,  98) = mat(k,  98) + lmat(k,  98)
         mat(k,  99) = lmat(k,  99)
         mat(k, 100) = lmat(k, 100)
         mat(k, 101) = lmat(k, 101)
         mat(k, 102) = mat(k, 102) + lmat(k, 102)
         mat(k, 105) = mat(k, 105) + lmat(k, 105)
         mat(k, 107) = mat(k, 107) + lmat(k, 107)
         mat(k, 108) = mat(k, 108) + lmat(k, 108)
         mat(k, 111) = mat(k, 111) + lmat(k, 111)
         mat(k, 112) = mat(k, 112) + lmat(k, 112)
         mat(k, 113) = mat(k, 113) + lmat(k, 113)
         mat(k, 115) = mat(k, 115) + lmat(k, 115)
         mat(k, 119) = mat(k, 119) + lmat(k, 119)
         mat(k, 120) = mat(k, 120) + lmat(k, 120)
         mat(k, 122) = mat(k, 122) + lmat(k, 122)
         mat(k, 123) = mat(k, 123) + lmat(k, 123)
         mat(k, 124) = mat(k, 124) + lmat(k, 124)
         mat(k, 126) = mat(k, 126) + lmat(k, 126)
         mat(k, 127) = mat(k, 127) + lmat(k, 127)
         mat(k, 128) = mat(k, 128) + lmat(k, 128)
         mat(k, 130) = mat(k, 130) + lmat(k, 130)
         mat(k, 131) = mat(k, 131) + lmat(k, 131)
         mat(k, 133) = mat(k, 133) + lmat(k, 133)
         mat(k, 135) = lmat(k, 135)
         mat(k, 136) = lmat(k, 136)
         mat(k, 137) = lmat(k, 137)
         mat(k, 138) = lmat(k, 138)
         mat(k, 139) = lmat(k, 139)
         mat(k, 140) = lmat(k, 140)
         mat(k, 141) = mat(k, 141) + lmat(k, 141)
         mat(k, 144) = lmat(k, 144)
         mat(k, 145) = lmat(k, 145)
         mat(k, 146) = lmat(k, 146)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 149) = mat(k, 149) + lmat(k, 149)
         mat(k, 151) = mat(k, 151) + lmat(k, 151)
         mat(k, 152) = mat(k, 152) + lmat(k, 152)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 154) = mat(k, 154) + lmat(k, 154)
         mat(k, 155) = mat(k, 155) + lmat(k, 155)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 159) = mat(k, 159) + lmat(k, 159)
         mat(k, 161) = mat(k, 161) + lmat(k, 161)
         mat(k, 162) = mat(k, 162) + lmat(k, 162)
         mat(k, 163) = mat(k, 163) + lmat(k, 163)
         mat(k, 164) = mat(k, 164) + lmat(k, 164)
         mat(k, 166) = mat(k, 166) + lmat(k, 166)
         mat(k, 168) = mat(k, 168) + lmat(k, 168)
         mat(k, 174) = lmat(k, 174)
         mat(k, 175) = lmat(k, 175)
         mat(k, 176) = lmat(k, 176)
         mat(k, 177) = lmat(k, 177)
         mat(k, 178) = mat(k, 178) + lmat(k, 178)
         mat(k, 183) = lmat(k, 183)
         mat(k, 184) = lmat(k, 184)
         mat(k, 185) = mat(k, 185) + lmat(k, 185)
         mat(k, 189) = mat(k, 189) + lmat(k, 189)
         mat(k, 190) = mat(k, 190) + lmat(k, 190)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 200) = mat(k, 200) + lmat(k, 200)
         mat(k, 205) = lmat(k, 205)
         mat(k, 206) = lmat(k, 206)
         mat(k, 207) = lmat(k, 207)
         mat(k, 208) = lmat(k, 208)
         mat(k, 209) = mat(k, 209) + lmat(k, 209)
         mat(k, 211) = mat(k, 211) + lmat(k, 211)
         mat(k, 214) = mat(k, 214) + lmat(k, 214)
         mat(k, 222) = mat(k, 222) + lmat(k, 222)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 228) = mat(k, 228) + lmat(k, 228)
         mat(k, 231) = mat(k, 231) + lmat(k, 231)
         mat(k, 232) = mat(k, 232) + lmat(k, 232)
         mat(k, 233) = mat(k, 233) + lmat(k, 233)
         mat(k, 235) = mat(k, 235) + lmat(k, 235)
         mat(k, 236) = mat(k, 236) + lmat(k, 236)
         mat(k, 237) = lmat(k, 237)
         mat(k, 239) = mat(k, 239) + lmat(k, 239)
         mat(k, 240) = lmat(k, 240)
         mat(k, 241) = lmat(k, 241)
         mat(k, 242) = lmat(k, 242)
         mat(k, 243) = lmat(k, 243)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 247) = lmat(k, 247)
         mat(k, 248) = lmat(k, 248)
         mat(k, 249) = lmat(k, 249)
         mat(k, 250) = mat(k, 250) + lmat(k, 250)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 254) = mat(k, 254) + lmat(k, 254)
         mat(k, 256) = mat(k, 256) + lmat(k, 256)
         mat(k, 257) = mat(k, 257) + lmat(k, 257)
         mat(k, 261) = mat(k, 261) + lmat(k, 261)
         mat(k, 262) = mat(k, 262) + lmat(k, 262)
         mat(k, 267) = mat(k, 267) + lmat(k, 267)
         mat(k, 268) = lmat(k, 268)
         mat(k, 269) = lmat(k, 269)
         mat(k, 271) = mat(k, 271) + lmat(k, 271)
         mat(k, 272) = lmat(k, 272)
         mat(k, 274) = mat(k, 274) + lmat(k, 274)
         mat(k, 275) = lmat(k, 275)
         mat(k, 276) = lmat(k, 276)
         mat(k, 277) = lmat(k, 277)
         mat(k, 278) = lmat(k, 278)
         mat(k, 279) = lmat(k, 279)
         mat(k, 280) = lmat(k, 280)
         mat(k, 281) = lmat(k, 281)
         mat(k, 282) = lmat(k, 282)
         mat(k, 283) = lmat(k, 283)
         mat(k, 284) = lmat(k, 284)
         mat(k, 285) = mat(k, 285) + lmat(k, 285)
         mat(k, 291) = lmat(k, 291)
         mat(k, 292) = lmat(k, 292)
         mat(k, 293) = lmat(k, 293)
         mat(k, 294) = mat(k, 294) + lmat(k, 294)
         mat(k, 300) = lmat(k, 300)
         mat(k, 301) = lmat(k, 301)
         mat(k, 302) = lmat(k, 302)
         mat(k, 303) = mat(k, 303) + lmat(k, 303)
         mat(k, 309) = mat(k, 309) + lmat(k, 309)
         mat(k, 310) = mat(k, 310) + lmat(k, 310)
         mat(k, 314) = mat(k, 314) + lmat(k, 314)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 321) = mat(k, 321) + lmat(k, 321)
         mat(k, 323) = lmat(k, 323)
         mat(k, 324) = lmat(k, 324)
         mat(k, 325) = mat(k, 325) + lmat(k, 325)
         mat(k, 326) = lmat(k, 326)
         mat(k, 327) = mat(k, 327) + lmat(k, 327)
         mat(k, 330) = mat(k, 330) + lmat(k, 330)
         mat(k, 331) = lmat(k, 331)
         mat(k, 332) = mat(k, 332) + lmat(k, 332)
         mat(k, 335) = mat(k, 335) + lmat(k, 335)
         mat(k, 338) = mat(k, 338) + lmat(k, 338)
         mat(k, 339) = lmat(k, 339)
         mat(k, 340) = mat(k, 340) + lmat(k, 340)
         mat(k, 342) = lmat(k, 342)
         mat(k, 343) = mat(k, 343) + lmat(k, 343)
         mat(k, 344) = lmat(k, 344)
         mat(k, 345) = mat(k, 345) + lmat(k, 345)
         mat(k, 347) = mat(k, 347) + lmat(k, 347)
         mat(k, 348) = mat(k, 348) + lmat(k, 348)
         mat(k, 349) = lmat(k, 349)
         mat(k, 350) = lmat(k, 350)
         mat(k, 351) = lmat(k, 351)
         mat(k, 352) = lmat(k, 352)
         mat(k, 353) = lmat(k, 353)
         mat(k, 354) = lmat(k, 354)
         mat(k, 355) = mat(k, 355) + lmat(k, 355)
         mat(k, 363) = mat(k, 363) + lmat(k, 363)
         mat(k, 364) = lmat(k, 364)
         mat(k, 366) = mat(k, 366) + lmat(k, 366)
         mat(k, 371) = mat(k, 371) + lmat(k, 371)
         mat(k, 374) = lmat(k, 374)
         mat(k, 376) = mat(k, 376) + lmat(k, 376)
         mat(k, 379) = mat(k, 379) + lmat(k, 379)
         mat(k, 382) = lmat(k, 382)
         mat(k, 385) = mat(k, 385) + lmat(k, 385)
         mat(k, 386) = lmat(k, 386)
         mat(k, 388) = lmat(k, 388)
         mat(k, 389) = mat(k, 389) + lmat(k, 389)
         mat(k, 390) = lmat(k, 390)
         mat(k, 391) = mat(k, 391) + lmat(k, 391)
         mat(k, 392) = lmat(k, 392)
         mat(k, 395) = mat(k, 395) + lmat(k, 395)
         mat(k, 396) = lmat(k, 396)
         mat(k, 397) = mat(k, 397) + lmat(k, 397)
         mat(k, 398) = lmat(k, 398)
         mat(k, 400) = lmat(k, 400)
         mat(k, 401) = mat(k, 401) + lmat(k, 401)
         mat(k, 402) = lmat(k, 402)
         mat(k, 403) = mat(k, 403) + lmat(k, 403)
         mat(k, 405) = lmat(k, 405)
         mat(k, 406) = mat(k, 406) + lmat(k, 406)
         mat(k, 407) = lmat(k, 407)
         mat(k, 409) = mat(k, 409) + lmat(k, 409)
         mat(k, 411) = lmat(k, 411)
         mat(k, 412) = mat(k, 412) + lmat(k, 412)
         mat(k, 413) = lmat(k, 413)
         mat(k, 415) = lmat(k, 415)
         mat(k, 416) = lmat(k, 416)
         mat(k, 417) = lmat(k, 417)
         mat(k, 418) = lmat(k, 418)
         mat(k, 419) = lmat(k, 419)
         mat(k, 420) = lmat(k, 420)
         mat(k, 421) = mat(k, 421) + lmat(k, 421)
         mat(k, 422) = mat(k, 422) + lmat(k, 422)
         mat(k, 423) = lmat(k, 423)
         mat(k, 424) = lmat(k, 424)
         mat(k, 425) = mat(k, 425) + lmat(k, 425)
         mat(k, 427) = mat(k, 427) + lmat(k, 427)
         mat(k, 429) = lmat(k, 429)
         mat(k, 430) = lmat(k, 430)
         mat(k, 431) = mat(k, 431) + lmat(k, 431)
         mat(k, 432) = lmat(k, 432)
         mat(k, 435) = mat(k, 435) + lmat(k, 435)
         mat(k, 441) = mat(k, 441) + lmat(k, 441)
         mat(k, 444) = mat(k, 444) + lmat(k, 444)
         mat(k, 445) = lmat(k, 445)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 448) = lmat(k, 448)
         mat(k, 449) = mat(k, 449) + lmat(k, 449)
         mat(k, 451) = mat(k, 451) + lmat(k, 451)
         mat(k, 453) = mat(k, 453) + lmat(k, 453)
         mat(k, 454) = mat(k, 454) + lmat(k, 454)
         mat(k, 460) = mat(k, 460) + lmat(k, 460)
         mat(k, 466) = lmat(k, 466)
         mat(k, 467) = lmat(k, 467)
         mat(k, 468) = lmat(k, 468)
         mat(k, 469) = mat(k, 469) + lmat(k, 469)
         mat(k, 473) = mat(k, 473) + lmat(k, 473)
         mat(k, 474) = mat(k, 474) + lmat(k, 474)
         mat(k, 481) = mat(k, 481) + lmat(k, 481)
         mat(k, 485) = mat(k, 485) + lmat(k, 485)
         mat(k, 487) = mat(k, 487) + lmat(k, 487)
         mat(k, 488) = mat(k, 488) + lmat(k, 488)
         mat(k, 489) = lmat(k, 489)
         mat(k, 490) = lmat(k, 490)
         mat(k, 492) = mat(k, 492) + lmat(k, 492)
         mat(k, 493) = lmat(k, 493)
         mat(k, 494) = lmat(k, 494)
         mat(k, 497) = mat(k, 497) + lmat(k, 497)
         mat(k, 505) = mat(k, 505) + lmat(k, 505)
         mat(k, 512) = mat(k, 512) + lmat(k, 512)
         mat(k, 514) = lmat(k, 514)
         mat(k, 515) = lmat(k, 515)
         mat(k, 517) = mat(k, 517) + lmat(k, 517)
         mat(k, 518) = mat(k, 518) + lmat(k, 518)
         mat(k, 522) = lmat(k, 522)
         mat(k, 524) = mat(k, 524) + lmat(k, 524)
         mat(k, 530) = mat(k, 530) + lmat(k, 530)
         mat(k, 531) = lmat(k, 531)
         mat(k, 532) = lmat(k, 532)
         mat(k, 534) = mat(k, 534) + lmat(k, 534)
         mat(k, 535) = lmat(k, 535)
         mat(k, 537) = lmat(k, 537)
         mat(k, 538) = mat(k, 538) + lmat(k, 538)
         mat(k, 540) = lmat(k, 540)
         mat(k, 542) = mat(k, 542) + lmat(k, 542)
         mat(k, 543) = lmat(k, 543)
         mat(k, 544) = lmat(k, 544)
         mat(k, 546) = mat(k, 546) + lmat(k, 546)
         mat(k, 551) = lmat(k, 551)
         mat(k, 554) = mat(k, 554) + lmat(k, 554)
         mat(k, 562) = mat(k, 562) + lmat(k, 562)
         mat(k, 563) = lmat(k, 563)
         mat(k, 564) = lmat(k, 564)
         mat(k, 565) = lmat(k, 565)
         mat(k, 566) = mat(k, 566) + lmat(k, 566)
         mat(k, 567) = mat(k, 567) + lmat(k, 567)
         mat(k, 568) = lmat(k, 568)
         mat(k, 570) = mat(k, 570) + lmat(k, 570)
         mat(k, 578) = mat(k, 578) + lmat(k, 578)
         mat(k, 582) = lmat(k, 582)
         mat(k, 587) = lmat(k, 587)
         mat(k, 588) = lmat(k, 588)
         mat(k, 589) = lmat(k, 589)
         mat(k, 590) = lmat(k, 590)
         mat(k, 591) = mat(k, 591) + lmat(k, 591)
         mat(k, 592) = mat(k, 592) + lmat(k, 592)
         mat(k, 595) = lmat(k, 595)
         mat(k, 600) = mat(k, 600) + lmat(k, 600)
         mat(k, 601) = mat(k, 601) + lmat(k, 601)
         mat(k, 606) = mat(k, 606) + lmat(k, 606)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 610) = lmat(k, 610)
         mat(k, 615) = lmat(k, 615)
         mat(k, 616) = mat(k, 616) + lmat(k, 616)
         mat(k, 617) = mat(k, 617) + lmat(k, 617)
         mat(k, 619) = mat(k, 619) + lmat(k, 619)
         mat(k, 620) = lmat(k, 620)
         mat(k, 623) = mat(k, 623) + lmat(k, 623)
         mat(k, 624) = mat(k, 624) + lmat(k, 624)
         mat(k, 625) = lmat(k, 625)
         mat(k, 627) = lmat(k, 627)
         mat(k, 628) = mat(k, 628) + lmat(k, 628)
         mat(k, 631) = mat(k, 631) + lmat(k, 631)
         mat(k, 637) = lmat(k, 637)
         mat(k, 638) = mat(k, 638) + lmat(k, 638)
         mat(k, 641) = mat(k, 641) + lmat(k, 641)
         mat(k, 642) = mat(k, 642) + lmat(k, 642)
         mat(k, 645) = lmat(k, 645)
         mat(k, 646) = mat(k, 646) + lmat(k, 646)
         mat(k, 647) = mat(k, 647) + lmat(k, 647)
         mat(k, 648) = mat(k, 648) + lmat(k, 648)
         mat(k, 649) = lmat(k, 649)
         mat(k, 650) = lmat(k, 650)
         mat(k, 651) = lmat(k, 651)
         mat(k, 652) = lmat(k, 652)
         mat(k, 653) = lmat(k, 653)
         mat(k, 655) = mat(k, 655) + lmat(k, 655)
         mat(k, 656) = lmat(k, 656)
         mat(k, 657) = lmat(k, 657)
         mat(k, 658) = lmat(k, 658)
         mat(k, 659) = lmat(k, 659)
         mat(k, 660) = lmat(k, 660)
         mat(k, 661) = mat(k, 661) + lmat(k, 661)
         mat(k, 666) = lmat(k, 666)
         mat(k, 668) = lmat(k, 668)
         mat(k, 669) = mat(k, 669) + lmat(k, 669)
         mat(k, 670) = lmat(k, 670)
         mat(k, 671) = lmat(k, 671)
         mat(k, 672) = mat(k, 672) + lmat(k, 672)
         mat(k, 676) = lmat(k, 676)
         mat(k, 677) = lmat(k, 677)
         mat(k, 679) = mat(k, 679) + lmat(k, 679)
         mat(k, 680) = lmat(k, 680)
         mat(k, 681) = lmat(k, 681)
         mat(k, 683) = mat(k, 683) + lmat(k, 683)
         mat(k, 690) = mat(k, 690) + lmat(k, 690)
         mat(k, 696) = mat(k, 696) + lmat(k, 696)
         mat(k, 698) = lmat(k, 698)
         mat(k, 699) = lmat(k, 699)
         mat(k, 700) = mat(k, 700) + lmat(k, 700)
         mat(k, 701) = lmat(k, 701)
         mat(k, 702) = lmat(k, 702)
         mat(k, 703) = lmat(k, 703)
         mat(k, 704) = lmat(k, 704)
         mat(k, 705) = lmat(k, 705)
         mat(k, 706) = lmat(k, 706)
         mat(k, 707) = mat(k, 707) + lmat(k, 707)
         mat(k, 712) = lmat(k, 712)
         mat(k, 714) = lmat(k, 714)
         mat(k, 716) = mat(k, 716) + lmat(k, 716)
         mat(k, 717) = lmat(k, 717)
         mat(k, 720) = mat(k, 720) + lmat(k, 720)
         mat(k, 727) = mat(k, 727) + lmat(k, 727)
         mat(k, 738) = mat(k, 738) + lmat(k, 738)
         mat(k, 754) = mat(k, 754) + lmat(k, 754)
         mat(k, 765) = mat(k, 765) + lmat(k, 765)
         mat(k, 774) = mat(k, 774) + lmat(k, 774)
         mat(k, 784) = mat(k, 784) + lmat(k, 784)
         mat(k, 785) = lmat(k, 785)
         mat(k, 787) = mat(k, 787) + lmat(k, 787)
         mat(k, 793) = mat(k, 793) + lmat(k, 793)
         mat(k, 803) = mat(k, 803) + lmat(k, 803)
         mat(k, 808) = mat(k, 808) + lmat(k, 808)
         mat(k, 809) = mat(k, 809) + lmat(k, 809)
         mat(k, 811) = mat(k, 811) + lmat(k, 811)
         mat(k, 815) = mat(k, 815) + lmat(k, 815)
         mat(k, 816) = mat(k, 816) + lmat(k, 816)
         mat(k, 817) = mat(k, 817) + lmat(k, 817)
         mat(k, 826) = mat(k, 826) + lmat(k, 826)
         mat(k, 834) = lmat(k, 834)
         mat(k, 835) = lmat(k, 835)
         mat(k, 836) = lmat(k, 836)
         mat(k, 837) = mat(k, 837) + lmat(k, 837)
         mat(k, 839) = lmat(k, 839)
         mat(k, 841) = lmat(k, 841)
         mat(k, 842) = mat(k, 842) + lmat(k, 842)
         mat(k, 844) = mat(k, 844) + lmat(k, 844)
         mat(k, 845) = lmat(k, 845)
         mat(k, 848) = lmat(k, 848)
         mat(k, 853) = mat(k, 853) + lmat(k, 853)
         mat(k, 869) = mat(k, 869) + lmat(k, 869)
         mat(k, 879) = mat(k, 879) + lmat(k, 879)
         mat(k, 897) = mat(k, 897) + lmat(k, 897)
         mat(k, 921) = mat(k, 921) + lmat(k, 921)
         mat(k, 932) = mat(k, 932) + lmat(k, 932)
         mat(k, 933) = mat(k, 933) + lmat(k, 933)
         mat(k, 934) = mat(k, 934) + lmat(k, 934)
         mat(k, 936) = mat(k, 936) + lmat(k, 936)
         mat(k, 937) = mat(k, 937) + lmat(k, 937)
         mat(k, 938) = mat(k, 938) + lmat(k, 938)
         mat(k, 940) = lmat(k, 940)
         mat(k, 944) = mat(k, 944) + lmat(k, 944)
         mat(k, 950) = mat(k, 950) + lmat(k, 950)
         mat(k, 952) = lmat(k, 952)
         mat(k, 954) = lmat(k, 954)
         mat(k, 956) = lmat(k, 956)
         mat(k, 958) = mat(k, 958) + lmat(k, 958)
         mat(k, 959) = mat(k, 959) + lmat(k, 959)
         mat(k, 961) = mat(k, 961) + lmat(k, 961)
         mat(k, 976) = mat(k, 976) + lmat(k, 976)
         mat(k, 996) = mat(k, 996) + lmat(k, 996)
         mat(k, 998) = lmat(k, 998)
         mat(k, 999) = lmat(k, 999)
         mat(k,1003) = lmat(k,1003)
         mat(k,1009) = mat(k,1009) + lmat(k,1009)
         mat(k,1026) = lmat(k,1026)
         mat(k,1030) = mat(k,1030) + lmat(k,1030)
         mat(k,1034) = lmat(k,1034)
         mat(k,1038) = mat(k,1038) + lmat(k,1038)
         mat(k,1040) = lmat(k,1040)
         mat(k,1042) = mat(k,1042) + lmat(k,1042)
         mat(k,1043) = lmat(k,1043)
         mat(k,1048) = lmat(k,1048)
         mat(k,1049) = lmat(k,1049)
         mat(k,1053) = mat(k,1053) + lmat(k,1053)
         mat(k,1054) = lmat(k,1054)
         mat(k,1056) = mat(k,1056) + lmat(k,1056)
         mat(k,1057) = mat(k,1057) + lmat(k,1057)
         mat(k,1063) = mat(k,1063) + lmat(k,1063)
         mat(k,1075) = mat(k,1075) + lmat(k,1075)
         mat(k,1076) = mat(k,1076) + lmat(k,1076)
         mat(k,1077) = mat(k,1077) + lmat(k,1077)
         mat(k,1078) = mat(k,1078) + lmat(k,1078)
         mat(k,1079) = mat(k,1079) + lmat(k,1079)
         mat(k,1080) = mat(k,1080) + lmat(k,1080)
         mat(k,1083) = mat(k,1083) + lmat(k,1083)
         mat(k,1084) = mat(k,1084) + lmat(k,1084)
         mat(k,1088) = mat(k,1088) + lmat(k,1088)
         mat(k,1093) = lmat(k,1093)
         mat(k,1094) = lmat(k,1094)
         mat(k,1095) = lmat(k,1095)
         mat(k,1096) = lmat(k,1096)
         mat(k,1097) = mat(k,1097) + lmat(k,1097)
         mat(k,1098) = lmat(k,1098)
         mat(k,1100) = lmat(k,1100)
         mat(k,1101) = lmat(k,1101)
         mat(k,1105) = mat(k,1105) + lmat(k,1105)
         mat(k,1106) = lmat(k,1106)
         mat(k,1107) = lmat(k,1107)
         mat(k,1110) = mat(k,1110) + lmat(k,1110)
         mat(k,1112) = lmat(k,1112)
         mat(k,1114) = mat(k,1114) + lmat(k,1114)
         mat(k,1115) = lmat(k,1115)
         mat(k,1124) = mat(k,1124) + lmat(k,1124)
         mat(k,1145) = mat(k,1145) + lmat(k,1145)
         mat(k,1161) = mat(k,1161) + lmat(k,1161)
         mat(k,1178) = mat(k,1178) + lmat(k,1178)
         mat(k,1198) = mat(k,1198) + lmat(k,1198)
         mat(k,1213) = mat(k,1213) + lmat(k,1213)
         mat(k,1214) = mat(k,1214) + lmat(k,1214)
         mat(k,1217) = mat(k,1217) + lmat(k,1217)
         mat(k,1218) = mat(k,1218) + lmat(k,1218)
         mat(k,1220) = mat(k,1220) + lmat(k,1220)
         mat(k,1221) = mat(k,1221) + lmat(k,1221)
         mat(k,1225) = mat(k,1225) + lmat(k,1225)
         mat(k,1226) = mat(k,1226) + lmat(k,1226)
         mat(k,1227) = mat(k,1227) + lmat(k,1227)
         mat(k,1231) = lmat(k,1231)
         mat(k,1235) = lmat(k,1235)
         mat(k,1236) = mat(k,1236) + lmat(k,1236)
         mat(k,1237) = mat(k,1237) + lmat(k,1237)
         mat(k,1248) = lmat(k,1248)
         mat(k,1262) = mat(k,1262) + lmat(k,1262)
         mat(k,1278) = lmat(k,1278)
         mat(k,1295) = mat(k,1295) + lmat(k,1295)
         mat(k,1305) = mat(k,1305) + lmat(k,1305)
         mat(k,1319) = mat(k,1319) + lmat(k,1319)
         mat(k,1334) = lmat(k,1334)
         mat(k,1336) = mat(k,1336) + lmat(k,1336)
         mat(k,1340) = mat(k,1340) + lmat(k,1340)
         mat(k,1342) = mat(k,1342) + lmat(k,1342)
         mat(k,1348) = lmat(k,1348)
         mat(k,1362) = mat(k,1362) + lmat(k,1362)
         mat(k,1394) = mat(k,1394) + lmat(k,1394)
         mat(k,1409) = mat(k,1409) + lmat(k,1409)
         mat(k,1423) = mat(k,1423) + lmat(k,1423)
         mat(k,1434) = lmat(k,1434)
         mat(k,1436) = lmat(k,1436)
         mat(k,1437) = mat(k,1437) + lmat(k,1437)
         mat(k,1438) = mat(k,1438) + lmat(k,1438)
         mat(k,1439) = mat(k,1439) + lmat(k,1439)
         mat(k,1442) = mat(k,1442) + lmat(k,1442)
         mat(k,1444) = mat(k,1444) + lmat(k,1444)
         mat(k,1446) = mat(k,1446) + lmat(k,1446)
         mat(k,1449) = lmat(k,1449)
         mat(k,1450) = mat(k,1450) + lmat(k,1450)
         mat(k,1453) = mat(k,1453) + lmat(k,1453)
         mat(k,1459) = mat(k,1459) + lmat(k,1459)
         mat(k,1471) = mat(k,1471) + lmat(k,1471)
         mat(k,1473) = mat(k,1473) + lmat(k,1473)
         mat(k,1481) = lmat(k,1481)
         mat(k,1490) = mat(k,1490) + lmat(k,1490)
         mat(k,1493) = mat(k,1493) + lmat(k,1493)
         mat(k,1506) = mat(k,1506) + lmat(k,1506)
         mat(k,1535) = mat(k,1535) + lmat(k,1535)
         mat(k,1547) = lmat(k,1547)
         mat(k,1632) = mat(k,1632) + lmat(k,1632)
         mat(k,1700) = mat(k,1700) + lmat(k,1700)
         mat(k,1712) = mat(k,1712) + lmat(k,1712)
         mat(k,1756) = mat(k,1756) + lmat(k,1756)
         mat(k,1760) = mat(k,1760) + lmat(k,1760)
         mat(k,1761) = mat(k,1761) + lmat(k,1761)
         mat(k,1769) = mat(k,1769) + lmat(k,1769)
         mat(k,1770) = mat(k,1770) + lmat(k,1770)
         mat(k,1795) = mat(k,1795) + lmat(k,1795)
         mat(k,1853) = mat(k,1853) + lmat(k,1853)
         mat(k,1862) = mat(k,1862) + lmat(k,1862)
         mat(k,1905) = mat(k,1905) + lmat(k,1905)
         mat(k,2013) = mat(k,2013) + lmat(k,2013)
         mat(k,2021) = mat(k,2021) + lmat(k,2021)
         mat(k,2024) = lmat(k,2024)
         mat(k,2025) = mat(k,2025) + lmat(k,2025)
         mat(k,2026) = lmat(k,2026)
         mat(k,2029) = mat(k,2029) + lmat(k,2029)
         mat(k,2038) = mat(k,2038) + lmat(k,2038)
         mat(k,2065) = mat(k,2065) + lmat(k,2065)
         mat(k,2066) = mat(k,2066) + lmat(k,2066)
         mat(k,2070) = mat(k,2070) + lmat(k,2070)
         mat(k,2105) = mat(k,2105) + lmat(k,2105)
         mat(k,2159) = mat(k,2159) + lmat(k,2159)
         mat(k,2168) = mat(k,2168) + lmat(k,2168)
         mat(k,2171) = mat(k,2171) + lmat(k,2171)
         mat(k,2179) = mat(k,2179) + lmat(k,2179)
         mat(k,2190) = mat(k,2190) + lmat(k,2190)
         mat(k,2192) = mat(k,2192) + lmat(k,2192)
         mat(k,2224) = mat(k,2224) + lmat(k,2224)
         mat(k,2227) = mat(k,2227) + lmat(k,2227)
         mat(k,2229) = mat(k,2229) + lmat(k,2229)
         mat(k,2237) = mat(k,2237) + lmat(k,2237)
         mat(k,2238) = mat(k,2238) + lmat(k,2238)
         mat(k,2266) = mat(k,2266) + lmat(k,2266)
         mat(k,2269) = mat(k,2269) + lmat(k,2269)
         mat(k,2277) = lmat(k,2277)
         mat(k,2280) = lmat(k,2280)
         mat(k,2283) = mat(k,2283) + lmat(k,2283)
         mat(k,2284) = mat(k,2284) + lmat(k,2284)
         mat(k,2295) = lmat(k,2295)
         mat(k,2296) = mat(k,2296) + lmat(k,2296)
         mat(k, 218) = 0._r8
         mat(k, 219) = 0._r8
         mat(k, 258) = 0._r8
         mat(k, 311) = 0._r8
         mat(k, 329) = 0._r8
         mat(k, 436) = 0._r8
         mat(k, 437) = 0._r8
         mat(k, 462) = 0._r8
         mat(k, 483) = 0._r8
         mat(k, 498) = 0._r8
         mat(k, 500) = 0._r8
         mat(k, 508) = 0._r8
         mat(k, 632) = 0._r8
         mat(k, 635) = 0._r8
         mat(k, 639) = 0._r8
         mat(k, 640) = 0._r8
         mat(k, 643) = 0._r8
         mat(k, 662) = 0._r8
         mat(k, 664) = 0._r8
         mat(k, 665) = 0._r8
         mat(k, 667) = 0._r8
         mat(k, 673) = 0._r8
         mat(k, 674) = 0._r8
         mat(k, 678) = 0._r8
         mat(k, 708) = 0._r8
         mat(k, 710) = 0._r8
         mat(k, 711) = 0._r8
         mat(k, 713) = 0._r8
         mat(k, 715) = 0._r8
         mat(k, 721) = 0._r8
         mat(k, 722) = 0._r8
         mat(k, 737) = 0._r8
         mat(k, 739) = 0._r8
         mat(k, 740) = 0._r8
         mat(k, 742) = 0._r8
         mat(k, 744) = 0._r8
         mat(k, 753) = 0._r8
         mat(k, 755) = 0._r8
         mat(k, 756) = 0._r8
         mat(k, 758) = 0._r8
         mat(k, 760) = 0._r8
         mat(k, 761) = 0._r8
         mat(k, 776) = 0._r8
         mat(k, 777) = 0._r8
         mat(k, 782) = 0._r8
         mat(k, 788) = 0._r8
         mat(k, 796) = 0._r8
         mat(k, 802) = 0._r8
         mat(k, 828) = 0._r8
         mat(k, 849) = 0._r8
         mat(k, 898) = 0._r8
         mat(k, 901) = 0._r8
         mat(k, 908) = 0._r8
         mat(k, 915) = 0._r8
         mat(k, 919) = 0._r8
         mat(k, 920) = 0._r8
         mat(k, 924) = 0._r8
         mat(k, 925) = 0._r8
         mat(k, 927) = 0._r8
         mat(k, 942) = 0._r8
         mat(k, 945) = 0._r8
         mat(k, 948) = 0._r8
         mat(k, 962) = 0._r8
         mat(k, 972) = 0._r8
         mat(k, 975) = 0._r8
         mat(k, 978) = 0._r8
         mat(k, 984) = 0._r8
         mat(k, 988) = 0._r8
         mat(k, 989) = 0._r8
         mat(k, 993) = 0._r8
         mat(k, 994) = 0._r8
         mat(k, 995) = 0._r8
         mat(k,1012) = 0._r8
         mat(k,1013) = 0._r8
         mat(k,1014) = 0._r8
         mat(k,1018) = 0._r8
         mat(k,1023) = 0._r8
         mat(k,1024) = 0._r8
         mat(k,1025) = 0._r8
         mat(k,1028) = 0._r8
         mat(k,1031) = 0._r8
         mat(k,1032) = 0._r8
         mat(k,1033) = 0._r8
         mat(k,1036) = 0._r8
         mat(k,1037) = 0._r8
         mat(k,1039) = 0._r8
         mat(k,1041) = 0._r8
         mat(k,1058) = 0._r8
         mat(k,1064) = 0._r8
         mat(k,1065) = 0._r8
         mat(k,1066) = 0._r8
         mat(k,1067) = 0._r8
         mat(k,1068) = 0._r8
         mat(k,1082) = 0._r8
         mat(k,1085) = 0._r8
         mat(k,1086) = 0._r8
         mat(k,1091) = 0._r8
         mat(k,1099) = 0._r8
         mat(k,1103) = 0._r8
         mat(k,1104) = 0._r8
         mat(k,1108) = 0._r8
         mat(k,1109) = 0._r8
         mat(k,1125) = 0._r8
         mat(k,1126) = 0._r8
         mat(k,1127) = 0._r8
         mat(k,1128) = 0._r8
         mat(k,1134) = 0._r8
         mat(k,1135) = 0._r8
         mat(k,1141) = 0._r8
         mat(k,1142) = 0._r8
         mat(k,1143) = 0._r8
         mat(k,1144) = 0._r8
         mat(k,1146) = 0._r8
         mat(k,1147) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1149) = 0._r8
         mat(k,1150) = 0._r8
         mat(k,1156) = 0._r8
         mat(k,1157) = 0._r8
         mat(k,1170) = 0._r8
         mat(k,1181) = 0._r8
         mat(k,1188) = 0._r8
         mat(k,1191) = 0._r8
         mat(k,1193) = 0._r8
         mat(k,1195) = 0._r8
         mat(k,1196) = 0._r8
         mat(k,1197) = 0._r8
         mat(k,1199) = 0._r8
         mat(k,1200) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1203) = 0._r8
         mat(k,1204) = 0._r8
         mat(k,1211) = 0._r8
         mat(k,1223) = 0._r8
         mat(k,1232) = 0._r8
         mat(k,1243) = 0._r8
         mat(k,1249) = 0._r8
         mat(k,1254) = 0._r8
         mat(k,1255) = 0._r8
         mat(k,1256) = 0._r8
         mat(k,1257) = 0._r8
         mat(k,1258) = 0._r8
         mat(k,1259) = 0._r8
         mat(k,1261) = 0._r8
         mat(k,1263) = 0._r8
         mat(k,1265) = 0._r8
         mat(k,1267) = 0._r8
         mat(k,1268) = 0._r8
         mat(k,1274) = 0._r8
         mat(k,1276) = 0._r8
         mat(k,1277) = 0._r8
         mat(k,1281) = 0._r8
         mat(k,1284) = 0._r8
         mat(k,1285) = 0._r8
         mat(k,1288) = 0._r8
         mat(k,1289) = 0._r8
         mat(k,1291) = 0._r8
         mat(k,1292) = 0._r8
         mat(k,1293) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1298) = 0._r8
         mat(k,1300) = 0._r8
         mat(k,1301) = 0._r8
         mat(k,1307) = 0._r8
         mat(k,1309) = 0._r8
         mat(k,1310) = 0._r8
         mat(k,1317) = 0._r8
         mat(k,1320) = 0._r8
         mat(k,1322) = 0._r8
         mat(k,1323) = 0._r8
         mat(k,1330) = 0._r8
         mat(k,1331) = 0._r8
         mat(k,1337) = 0._r8
         mat(k,1341) = 0._r8
         mat(k,1343) = 0._r8
         mat(k,1344) = 0._r8
         mat(k,1346) = 0._r8
         mat(k,1347) = 0._r8
         mat(k,1352) = 0._r8
         mat(k,1353) = 0._r8
         mat(k,1354) = 0._r8
         mat(k,1359) = 0._r8
         mat(k,1360) = 0._r8
         mat(k,1361) = 0._r8
         mat(k,1364) = 0._r8
         mat(k,1373) = 0._r8
         mat(k,1374) = 0._r8
         mat(k,1395) = 0._r8
         mat(k,1396) = 0._r8
         mat(k,1398) = 0._r8
         mat(k,1405) = 0._r8
         mat(k,1406) = 0._r8
         mat(k,1410) = 0._r8
         mat(k,1412) = 0._r8
         mat(k,1415) = 0._r8
         mat(k,1416) = 0._r8
         mat(k,1424) = 0._r8
         mat(k,1426) = 0._r8
         mat(k,1427) = 0._r8
         mat(k,1431) = 0._r8
         mat(k,1432) = 0._r8
         mat(k,1440) = 0._r8
         mat(k,1443) = 0._r8
         mat(k,1452) = 0._r8
         mat(k,1454) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1457) = 0._r8
         mat(k,1458) = 0._r8
         mat(k,1460) = 0._r8
         mat(k,1461) = 0._r8
         mat(k,1463) = 0._r8
         mat(k,1468) = 0._r8
         mat(k,1469) = 0._r8
         mat(k,1470) = 0._r8
         mat(k,1472) = 0._r8
         mat(k,1475) = 0._r8
         mat(k,1476) = 0._r8
         mat(k,1477) = 0._r8
         mat(k,1478) = 0._r8
         mat(k,1479) = 0._r8
         mat(k,1480) = 0._r8
         mat(k,1482) = 0._r8
         mat(k,1491) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1494) = 0._r8
         mat(k,1496) = 0._r8
         mat(k,1498) = 0._r8
         mat(k,1500) = 0._r8
         mat(k,1503) = 0._r8
         mat(k,1504) = 0._r8
         mat(k,1507) = 0._r8
         mat(k,1533) = 0._r8
         mat(k,1537) = 0._r8
         mat(k,1546) = 0._r8
         mat(k,1614) = 0._r8
         mat(k,1630) = 0._r8
         mat(k,1646) = 0._r8
         mat(k,1649) = 0._r8
         mat(k,1660) = 0._r8
         mat(k,1661) = 0._r8
         mat(k,1683) = 0._r8
         mat(k,1699) = 0._r8
         mat(k,1719) = 0._r8
         mat(k,1724) = 0._r8
         mat(k,1725) = 0._r8
         mat(k,1735) = 0._r8
         mat(k,1737) = 0._r8
         mat(k,1740) = 0._r8
         mat(k,1741) = 0._r8
         mat(k,1746) = 0._r8
         mat(k,1753) = 0._r8
         mat(k,1754) = 0._r8
         mat(k,1755) = 0._r8
         mat(k,1757) = 0._r8
         mat(k,1758) = 0._r8
         mat(k,1762) = 0._r8
         mat(k,1765) = 0._r8
         mat(k,1766) = 0._r8
         mat(k,1767) = 0._r8
         mat(k,1768) = 0._r8
         mat(k,1771) = 0._r8
         mat(k,1817) = 0._r8
         mat(k,1846) = 0._r8
         mat(k,1847) = 0._r8
         mat(k,1848) = 0._r8
         mat(k,1850) = 0._r8
         mat(k,1860) = 0._r8
         mat(k,1863) = 0._r8
         mat(k,1874) = 0._r8
         mat(k,1898) = 0._r8
         mat(k,1899) = 0._r8
         mat(k,1900) = 0._r8
         mat(k,1901) = 0._r8
         mat(k,1902) = 0._r8
         mat(k,1903) = 0._r8
         mat(k,1910) = 0._r8
         mat(k,1911) = 0._r8
         mat(k,1913) = 0._r8
         mat(k,1914) = 0._r8
         mat(k,1930) = 0._r8
         mat(k,1947) = 0._r8
         mat(k,1948) = 0._r8
         mat(k,1976) = 0._r8
         mat(k,1980) = 0._r8
         mat(k,1982) = 0._r8
         mat(k,1984) = 0._r8
         mat(k,1986) = 0._r8
         mat(k,1988) = 0._r8
         mat(k,1994) = 0._r8
         mat(k,1999) = 0._r8
         mat(k,2004) = 0._r8
         mat(k,2006) = 0._r8
         mat(k,2008) = 0._r8
         mat(k,2028) = 0._r8
         mat(k,2031) = 0._r8
         mat(k,2032) = 0._r8
         mat(k,2035) = 0._r8
         mat(k,2036) = 0._r8
         mat(k,2039) = 0._r8
         mat(k,2041) = 0._r8
         mat(k,2043) = 0._r8
         mat(k,2055) = 0._r8
         mat(k,2056) = 0._r8
         mat(k,2058) = 0._r8
         mat(k,2060) = 0._r8
         mat(k,2067) = 0._r8
         mat(k,2071) = 0._r8
         mat(k,2081) = 0._r8
         mat(k,2085) = 0._r8
         mat(k,2086) = 0._r8
         mat(k,2087) = 0._r8
         mat(k,2088) = 0._r8
         mat(k,2090) = 0._r8
         mat(k,2095) = 0._r8
         mat(k,2096) = 0._r8
         mat(k,2097) = 0._r8
         mat(k,2100) = 0._r8
         mat(k,2108) = 0._r8
         mat(k,2109) = 0._r8
         mat(k,2110) = 0._r8
         mat(k,2123) = 0._r8
         mat(k,2128) = 0._r8
         mat(k,2131) = 0._r8
         mat(k,2136) = 0._r8
         mat(k,2137) = 0._r8
         mat(k,2140) = 0._r8
         mat(k,2141) = 0._r8
         mat(k,2142) = 0._r8
         mat(k,2144) = 0._r8
         mat(k,2148) = 0._r8
         mat(k,2149) = 0._r8
         mat(k,2150) = 0._r8
         mat(k,2152) = 0._r8
         mat(k,2157) = 0._r8
         mat(k,2169) = 0._r8
         mat(k,2172) = 0._r8
         mat(k,2178) = 0._r8
         mat(k,2181) = 0._r8
         mat(k,2184) = 0._r8
         mat(k,2185) = 0._r8
         mat(k,2186) = 0._r8
         mat(k,2187) = 0._r8
         mat(k,2188) = 0._r8
         mat(k,2191) = 0._r8
         mat(k,2193) = 0._r8
         mat(k,2206) = 0._r8
         mat(k,2209) = 0._r8
         mat(k,2210) = 0._r8
         mat(k,2214) = 0._r8
         mat(k,2215) = 0._r8
         mat(k,2216) = 0._r8
         mat(k,2217) = 0._r8
         mat(k,2221) = 0._r8
         mat(k,2222) = 0._r8
         mat(k,2223) = 0._r8
         mat(k,2226) = 0._r8
         mat(k,2230) = 0._r8
         mat(k,2232) = 0._r8
         mat(k,2234) = 0._r8
         mat(k,2236) = 0._r8
         mat(k,2239) = 0._r8
         mat(k,2241) = 0._r8
         mat(k,2247) = 0._r8
         mat(k,2253) = 0._r8
         mat(k,2255) = 0._r8
         mat(k,2257) = 0._r8
         mat(k,2261) = 0._r8
         mat(k,2270) = 0._r8
         mat(k,2274) = 0._r8
         mat(k,2276) = 0._r8
         mat(k,2278) = 0._r8
         mat(k,2279) = 0._r8
         mat(k,2281) = 0._r8
         mat(k,2282) = 0._r8
         mat(k,2285) = 0._r8
         mat(k,2286) = 0._r8
         mat(k,2287) = 0._r8
         mat(k,2288) = 0._r8
         mat(k,2289) = 0._r8
         mat(k,2290) = 0._r8
         mat(k,2291) = 0._r8
         mat(k,2292) = 0._r8
         mat(k,2293) = 0._r8
         mat(k,2294) = 0._r8
         mat(k,   1) = mat(k,   1) - dti(k)
         mat(k,   2) = mat(k,   2) - dti(k)
         mat(k,   3) = mat(k,   3) - dti(k)
         mat(k,   4) = mat(k,   4) - dti(k)
         mat(k,   5) = mat(k,   5) - dti(k)
         mat(k,   6) = mat(k,   6) - dti(k)
         mat(k,   7) = mat(k,   7) - dti(k)
         mat(k,   8) = mat(k,   8) - dti(k)
         mat(k,   9) = mat(k,   9) - dti(k)
         mat(k,  10) = mat(k,  10) - dti(k)
         mat(k,  11) = mat(k,  11) - dti(k)
         mat(k,  12) = mat(k,  12) - dti(k)
         mat(k,  13) = mat(k,  13) - dti(k)
         mat(k,  14) = mat(k,  14) - dti(k)
         mat(k,  15) = mat(k,  15) - dti(k)
         mat(k,  16) = mat(k,  16) - dti(k)
         mat(k,  17) = mat(k,  17) - dti(k)
         mat(k,  18) = mat(k,  18) - dti(k)
         mat(k,  19) = mat(k,  19) - dti(k)
         mat(k,  20) = mat(k,  20) - dti(k)
         mat(k,  21) = mat(k,  21) - dti(k)
         mat(k,  22) = mat(k,  22) - dti(k)
         mat(k,  23) = mat(k,  23) - dti(k)
         mat(k,  24) = mat(k,  24) - dti(k)
         mat(k,  25) = mat(k,  25) - dti(k)
         mat(k,  26) = mat(k,  26) - dti(k)
         mat(k,  27) = mat(k,  27) - dti(k)
         mat(k,  28) = mat(k,  28) - dti(k)
         mat(k,  29) = mat(k,  29) - dti(k)
         mat(k,  30) = mat(k,  30) - dti(k)
         mat(k,  31) = mat(k,  31) - dti(k)
         mat(k,  32) = mat(k,  32) - dti(k)
         mat(k,  33) = mat(k,  33) - dti(k)
         mat(k,  34) = mat(k,  34) - dti(k)
         mat(k,  35) = mat(k,  35) - dti(k)
         mat(k,  36) = mat(k,  36) - dti(k)
         mat(k,  37) = mat(k,  37) - dti(k)
         mat(k,  38) = mat(k,  38) - dti(k)
         mat(k,  39) = mat(k,  39) - dti(k)
         mat(k,  40) = mat(k,  40) - dti(k)
         mat(k,  41) = mat(k,  41) - dti(k)
         mat(k,  42) = mat(k,  42) - dti(k)
         mat(k,  43) = mat(k,  43) - dti(k)
         mat(k,  44) = mat(k,  44) - dti(k)
         mat(k,  45) = mat(k,  45) - dti(k)
         mat(k,  46) = mat(k,  46) - dti(k)
         mat(k,  47) = mat(k,  47) - dti(k)
         mat(k,  48) = mat(k,  48) - dti(k)
         mat(k,  54) = mat(k,  54) - dti(k)
         mat(k,  60) = mat(k,  60) - dti(k)
         mat(k,  66) = mat(k,  66) - dti(k)
         mat(k,  72) = mat(k,  72) - dti(k)
         mat(k,  78) = mat(k,  78) - dti(k)
         mat(k,  80) = mat(k,  80) - dti(k)
         mat(k,  86) = mat(k,  86) - dti(k)
         mat(k,  92) = mat(k,  92) - dti(k)
         mat(k,  98) = mat(k,  98) - dti(k)
         mat(k,  99) = mat(k,  99) - dti(k)
         mat(k, 102) = mat(k, 102) - dti(k)
         mat(k, 105) = mat(k, 105) - dti(k)
         mat(k, 108) = mat(k, 108) - dti(k)
         mat(k, 111) = mat(k, 111) - dti(k)
         mat(k, 115) = mat(k, 115) - dti(k)
         mat(k, 119) = mat(k, 119) - dti(k)
         mat(k, 123) = mat(k, 123) - dti(k)
         mat(k, 127) = mat(k, 127) - dti(k)
         mat(k, 131) = mat(k, 131) - dti(k)
         mat(k, 135) = mat(k, 135) - dti(k)
         mat(k, 138) = mat(k, 138) - dti(k)
         mat(k, 141) = mat(k, 141) - dti(k)
         mat(k, 144) = mat(k, 144) - dti(k)
         mat(k, 147) = mat(k, 147) - dti(k)
         mat(k, 152) = mat(k, 152) - dti(k)
         mat(k, 157) = mat(k, 157) - dti(k)
         mat(k, 162) = mat(k, 162) - dti(k)
         mat(k, 168) = mat(k, 168) - dti(k)
         mat(k, 174) = mat(k, 174) - dti(k)
         mat(k, 178) = mat(k, 178) - dti(k)
         mat(k, 183) = mat(k, 183) - dti(k)
         mat(k, 185) = mat(k, 185) - dti(k)
         mat(k, 189) = mat(k, 189) - dti(k)
         mat(k, 193) = mat(k, 193) - dti(k)
         mat(k, 200) = mat(k, 200) - dti(k)
         mat(k, 205) = mat(k, 205) - dti(k)
         mat(k, 209) = mat(k, 209) - dti(k)
         mat(k, 214) = mat(k, 214) - dti(k)
         mat(k, 222) = mat(k, 222) - dti(k)
         mat(k, 227) = mat(k, 227) - dti(k)
         mat(k, 232) = mat(k, 232) - dti(k)
         mat(k, 236) = mat(k, 236) - dti(k)
         mat(k, 241) = mat(k, 241) - dti(k)
         mat(k, 244) = mat(k, 244) - dti(k)
         mat(k, 247) = mat(k, 247) - dti(k)
         mat(k, 250) = mat(k, 250) - dti(k)
         mat(k, 253) = mat(k, 253) - dti(k)
         mat(k, 257) = mat(k, 257) - dti(k)
         mat(k, 262) = mat(k, 262) - dti(k)
         mat(k, 267) = mat(k, 267) - dti(k)
         mat(k, 271) = mat(k, 271) - dti(k)
         mat(k, 275) = mat(k, 275) - dti(k)
         mat(k, 279) = mat(k, 279) - dti(k)
         mat(k, 285) = mat(k, 285) - dti(k)
         mat(k, 291) = mat(k, 291) - dti(k)
         mat(k, 294) = mat(k, 294) - dti(k)
         mat(k, 300) = mat(k, 300) - dti(k)
         mat(k, 303) = mat(k, 303) - dti(k)
         mat(k, 310) = mat(k, 310) - dti(k)
         mat(k, 316) = mat(k, 316) - dti(k)
         mat(k, 321) = mat(k, 321) - dti(k)
         mat(k, 327) = mat(k, 327) - dti(k)
         mat(k, 332) = mat(k, 332) - dti(k)
         mat(k, 335) = mat(k, 335) - dti(k)
         mat(k, 340) = mat(k, 340) - dti(k)
         mat(k, 345) = mat(k, 345) - dti(k)
         mat(k, 350) = mat(k, 350) - dti(k)
         mat(k, 355) = mat(k, 355) - dti(k)
         mat(k, 363) = mat(k, 363) - dti(k)
         mat(k, 371) = mat(k, 371) - dti(k)
         mat(k, 379) = mat(k, 379) - dti(k)
         mat(k, 385) = mat(k, 385) - dti(k)
         mat(k, 391) = mat(k, 391) - dti(k)
         mat(k, 397) = mat(k, 397) - dti(k)
         mat(k, 403) = mat(k, 403) - dti(k)
         mat(k, 409) = mat(k, 409) - dti(k)
         mat(k, 415) = mat(k, 415) - dti(k)
         mat(k, 421) = mat(k, 421) - dti(k)
         mat(k, 427) = mat(k, 427) - dti(k)
         mat(k, 435) = mat(k, 435) - dti(k)
         mat(k, 441) = mat(k, 441) - dti(k)
         mat(k, 447) = mat(k, 447) - dti(k)
         mat(k, 453) = mat(k, 453) - dti(k)
         mat(k, 460) = mat(k, 460) - dti(k)
         mat(k, 466) = mat(k, 466) - dti(k)
         mat(k, 469) = mat(k, 469) - dti(k)
         mat(k, 473) = mat(k, 473) - dti(k)
         mat(k, 481) = mat(k, 481) - dti(k)
         mat(k, 488) = mat(k, 488) - dti(k)
         mat(k, 497) = mat(k, 497) - dti(k)
         mat(k, 505) = mat(k, 505) - dti(k)
         mat(k, 512) = mat(k, 512) - dti(k)
         mat(k, 517) = mat(k, 517) - dti(k)
         mat(k, 524) = mat(k, 524) - dti(k)
         mat(k, 530) = mat(k, 530) - dti(k)
         mat(k, 538) = mat(k, 538) - dti(k)
         mat(k, 546) = mat(k, 546) - dti(k)
         mat(k, 554) = mat(k, 554) - dti(k)
         mat(k, 562) = mat(k, 562) - dti(k)
         mat(k, 570) = mat(k, 570) - dti(k)
         mat(k, 578) = mat(k, 578) - dti(k)
         mat(k, 587) = mat(k, 587) - dti(k)
         mat(k, 591) = mat(k, 591) - dti(k)
         mat(k, 600) = mat(k, 600) - dti(k)
         mat(k, 607) = mat(k, 607) - dti(k)
         mat(k, 616) = mat(k, 616) - dti(k)
         mat(k, 623) = mat(k, 623) - dti(k)
         mat(k, 631) = mat(k, 631) - dti(k)
         mat(k, 638) = mat(k, 638) - dti(k)
         mat(k, 648) = mat(k, 648) - dti(k)
         mat(k, 661) = mat(k, 661) - dti(k)
         mat(k, 672) = mat(k, 672) - dti(k)
         mat(k, 683) = mat(k, 683) - dti(k)
         mat(k, 690) = mat(k, 690) - dti(k)
         mat(k, 696) = mat(k, 696) - dti(k)
         mat(k, 707) = mat(k, 707) - dti(k)
         mat(k, 720) = mat(k, 720) - dti(k)
         mat(k, 727) = mat(k, 727) - dti(k)
         mat(k, 738) = mat(k, 738) - dti(k)
         mat(k, 754) = mat(k, 754) - dti(k)
         mat(k, 765) = mat(k, 765) - dti(k)
         mat(k, 774) = mat(k, 774) - dti(k)
         mat(k, 784) = mat(k, 784) - dti(k)
         mat(k, 793) = mat(k, 793) - dti(k)
         mat(k, 803) = mat(k, 803) - dti(k)
         mat(k, 808) = mat(k, 808) - dti(k)
         mat(k, 815) = mat(k, 815) - dti(k)
         mat(k, 826) = mat(k, 826) - dti(k)
         mat(k, 834) = mat(k, 834) - dti(k)
         mat(k, 837) = mat(k, 837) - dti(k)
         mat(k, 844) = mat(k, 844) - dti(k)
         mat(k, 853) = mat(k, 853) - dti(k)
         mat(k, 869) = mat(k, 869) - dti(k)
         mat(k, 879) = mat(k, 879) - dti(k)
         mat(k, 897) = mat(k, 897) - dti(k)
         mat(k, 921) = mat(k, 921) - dti(k)
         mat(k, 933) = mat(k, 933) - dti(k)
         mat(k, 944) = mat(k, 944) - dti(k)
         mat(k, 950) = mat(k, 950) - dti(k)
         mat(k, 958) = mat(k, 958) - dti(k)
         mat(k, 976) = mat(k, 976) - dti(k)
         mat(k, 996) = mat(k, 996) - dti(k)
         mat(k,1009) = mat(k,1009) - dti(k)
         mat(k,1030) = mat(k,1030) - dti(k)
         mat(k,1042) = mat(k,1042) - dti(k)
         mat(k,1053) = mat(k,1053) - dti(k)
         mat(k,1063) = mat(k,1063) - dti(k)
         mat(k,1077) = mat(k,1077) - dti(k)
         mat(k,1088) = mat(k,1088) - dti(k)
         mat(k,1097) = mat(k,1097) - dti(k)
         mat(k,1110) = mat(k,1110) - dti(k)
         mat(k,1124) = mat(k,1124) - dti(k)
         mat(k,1145) = mat(k,1145) - dti(k)
         mat(k,1161) = mat(k,1161) - dti(k)
         mat(k,1178) = mat(k,1178) - dti(k)
         mat(k,1198) = mat(k,1198) - dti(k)
         mat(k,1214) = mat(k,1214) - dti(k)
         mat(k,1226) = mat(k,1226) - dti(k)
         mat(k,1237) = mat(k,1237) - dti(k)
         mat(k,1262) = mat(k,1262) - dti(k)
         mat(k,1295) = mat(k,1295) - dti(k)
         mat(k,1319) = mat(k,1319) - dti(k)
         mat(k,1340) = mat(k,1340) - dti(k)
         mat(k,1362) = mat(k,1362) - dti(k)
         mat(k,1394) = mat(k,1394) - dti(k)
         mat(k,1409) = mat(k,1409) - dti(k)
         mat(k,1423) = mat(k,1423) - dti(k)
         mat(k,1438) = mat(k,1438) - dti(k)
         mat(k,1453) = mat(k,1453) - dti(k)
         mat(k,1471) = mat(k,1471) - dti(k)
         mat(k,1493) = mat(k,1493) - dti(k)
         mat(k,1535) = mat(k,1535) - dti(k)
         mat(k,1700) = mat(k,1700) - dti(k)
         mat(k,1760) = mat(k,1760) - dti(k)
         mat(k,1853) = mat(k,1853) - dti(k)
         mat(k,1905) = mat(k,1905) - dti(k)
         mat(k,2013) = mat(k,2013) - dti(k)
         mat(k,2038) = mat(k,2038) - dti(k)
         mat(k,2065) = mat(k,2065) - dti(k)
         mat(k,2105) = mat(k,2105) - dti(k)
         mat(k,2168) = mat(k,2168) - dti(k)
         mat(k,2192) = mat(k,2192) - dti(k)
         mat(k,2237) = mat(k,2237) - dti(k)
         mat(k,2269) = mat(k,2269) - dti(k)
         mat(k,2296) = mat(k,2296) - dti(k)
      end do

      end subroutine nlnmat_finit

      subroutine     nlnmat( avec_len, mat, y, rxt, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer, intent(in) ::  avec_len
      real(r8), intent(in)    ::  dti(veclen)
      real(r8), intent(in)    ::  lmat(veclen,nzcnt)
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)

      call     nlnmat01( avec_len, mat, y, rxt )
      call     nlnmat02( avec_len, mat, y, rxt )
      call     nlnmat03( avec_len, mat, y, rxt )
      call     nlnmat04( avec_len, mat, y, rxt )
      call     nlnmat05( avec_len, mat, y, rxt )
      call     nlnmat06( avec_len, mat, y, rxt )
      call     nlnmat07( avec_len, mat, y, rxt )
      call     nlnmat08( avec_len, mat, y, rxt )
      call     nlnmat09( avec_len, mat, y, rxt )
      call     nlnmat10( avec_len, mat, y, rxt )
      call     nlnmat_finit( avec_len, mat, lmat, dti )

      end subroutine nlnmat

      end module mo_nln_matrix

