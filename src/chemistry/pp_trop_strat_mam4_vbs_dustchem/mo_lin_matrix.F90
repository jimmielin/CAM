
      module mo_lin_matrix

      use chem_mods, only: veclen
      private
      public :: linmat

      contains

      subroutine linmat01( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
!       ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer,  intent(in)    ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(in)    ::  het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)
!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k


      do k = 1,avec_len

         mat(k,638) = -( rxt(k,19) + het_rates(k,1) )

         mat(k,672) = -( rxt(k,20) + het_rates(k,2) )

         mat(k,1) = -( het_rates(k,3) )

         mat(k,2) = -( het_rates(k,4) )

         mat(k,3) = -( het_rates(k,5) )

         mat(k,976) = -( het_rates(k,6) )

         mat(k,168) = -( het_rates(k,7) )

         mat(k,397) = -( rxt(k,21) + het_rates(k,8) )

         mat(k,174) = -( rxt(k,22) + het_rates(k,9) )

         mat(k,415) = -( rxt(k,23) + het_rates(k,10) )

         mat(k,466) = -( rxt(k,24) + het_rates(k,11) )
         mat(k,398) = .500_r8*rxt(k,21)
         mat(k,175) = rxt(k,22)
         mat(k,659) = .200_r8*rxt(k,70)
         mat(k,705) = .060_r8*rxt(k,72)

         mat(k,300) = -( rxt(k,25) + het_rates(k,12) )
         mat(k,658) = .200_r8*rxt(k,70)
         mat(k,703) = .200_r8*rxt(k,72)

         mat(k,587) = -( rxt(k,26) + het_rates(k,13) )
         mat(k,237) = rxt(k,46)
         mat(k,1026) = rxt(k,56)
         mat(k,660) = .200_r8*rxt(k,70)
         mat(k,706) = .150_r8*rxt(k,72)

         mat(k,350) = -( rxt(k,27) + het_rates(k,14) )
         mat(k,704) = .210_r8*rxt(k,72)

         mat(k,250) = -( het_rates(k,15) )

         mat(k,355) = -( het_rates(k,16) )

         mat(k,1423) = -( het_rates(k,17) )
         mat(k,242) = rxt(k,74)
         mat(k,1490) = rxt(k,75)
         mat(k,564) = rxt(k,77)
         mat(k,149) = rxt(k,79)
         mat(k,155) = rxt(k,80)
         mat(k,474) = 2.000_r8*rxt(k,86)
         mat(k,592) = rxt(k,87)
         mat(k,454) = 3.000_r8*rxt(k,90)
         mat(k,113) = 2.000_r8*rxt(k,98)
         mat(k,816) = rxt(k,99)
         mat(k,785) = rxt(k,105)

         mat(k,241) = -( rxt(k,74) + het_rates(k,18) )

         mat(k,1493) = -( rxt(k,75) + het_rates(k,19) )
         mat(k,566) = rxt(k,76)

         mat(k,562) = -( rxt(k,76) + rxt(k,77) + rxt(k,531) + rxt(k,534) + rxt(k,539) &
                 + het_rates(k,20) )

         mat(k,4) = -( het_rates(k,21) )

         mat(k,244) = -( het_rates(k,22) )
         mat(k,326) = rxt(k,28)

         mat(k,327) = -( rxt(k,28) + het_rates(k,23) )

         mat(k,285) = -( het_rates(k,24) )

         mat(k,554) = -( het_rates(k,25) )

         mat(k,262) = -( het_rates(k,26) )

         mat(k,345) = -( rxt(k,29) + het_rates(k,27) )

         mat(k,294) = -( het_rates(k,28) )

         mat(k,1009) = -( het_rates(k,29) )
         mat(k,1334) = .700_r8*rxt(k,55)

         mat(k,409) = -( rxt(k,30) + het_rates(k,30) )

         mat(k,115) = -( het_rates(k,31) )

         mat(k,271) = -( rxt(k,31) + het_rates(k,32) )

         mat(k,105) = -( rxt(k,78) + het_rates(k,33) )

         mat(k,147) = -( rxt(k,79) + het_rates(k,34) )

         mat(k,152) = -( rxt(k,80) + het_rates(k,35) )

         mat(k,119) = -( rxt(k,81) + het_rates(k,36) )

         mat(k,157) = -( rxt(k,82) + het_rates(k,37) )

         mat(k,123) = -( rxt(k,83) + het_rates(k,38) )

         mat(k,162) = -( rxt(k,84) + het_rates(k,39) )

         mat(k,127) = -( rxt(k,85) + het_rates(k,40) )

         mat(k,473) = -( rxt(k,86) + het_rates(k,41) )

         mat(k,2038) = -( rxt(k,32) + rxt(k,33) + rxt(k,500) + het_rates(k,42) )
         mat(k,646) = .100_r8*rxt(k,19)
         mat(k,681) = .100_r8*rxt(k,20)
         mat(k,451) = rxt(k,38)
         mat(k,1446) = .180_r8*rxt(k,39)
         mat(k,1057) = rxt(k,43)
         mat(k,1106) = .330_r8*rxt(k,45)
         mat(k,1115) = rxt(k,47)
         mat(k,702) = rxt(k,49)
         mat(k,1221) = 1.340_r8*rxt(k,50)
         mat(k,841) = rxt(k,57)
         mat(k,544) = rxt(k,62)
         mat(k,407) = rxt(k,63)
         mat(k,657) = .375_r8*rxt(k,65)
         mat(k,494) = .400_r8*rxt(k,67)
         mat(k,1084) = .680_r8*rxt(k,69)
         mat(k,445) = rxt(k,269)
         mat(k,278) = 2.000_r8*rxt(k,299)

         mat(k,591) = -( rxt(k,87) + het_rates(k,43) )

         mat(k,131) = -( rxt(k,88) + het_rates(k,44) )

         mat(k,1042) = -( rxt(k,34) + het_rates(k,45) )
         mat(k,642) = .400_r8*rxt(k,19)
         mat(k,677) = .400_r8*rxt(k,20)
         mat(k,347) = rxt(k,29)
         mat(k,1094) = .330_r8*rxt(k,45)
         mat(k,323) = rxt(k,53)
         mat(k,540) = rxt(k,62)

         mat(k,371) = -( rxt(k,89) + het_rates(k,46) )

         mat(k,108) = -( het_rates(k,47) )

         mat(k,950) = -( rxt(k,35) + het_rates(k,48) )
         mat(k,641) = .250_r8*rxt(k,19)
         mat(k,676) = .250_r8*rxt(k,20)
         mat(k,411) = .820_r8*rxt(k,30)
         mat(k,1093) = .170_r8*rxt(k,45)
         mat(k,650) = .300_r8*rxt(k,65)
         mat(k,489) = .050_r8*rxt(k,67)
         mat(k,1076) = .500_r8*rxt(k,69)

         mat(k,1226) = -( rxt(k,36) + het_rates(k,49) )
         mat(k,418) = .180_r8*rxt(k,23)
         mat(k,352) = rxt(k,27)
         mat(k,668) = .400_r8*rxt(k,70)
         mat(k,714) = .540_r8*rxt(k,72)
         mat(k,430) = .510_r8*rxt(k,73)

         mat(k,690) = -( het_rates(k,50) )

         mat(k,616) = -( rxt(k,37) + het_rates(k,51) )

         mat(k,803) = -( het_rates(k,52) )

         mat(k,447) = -( rxt(k,38) + het_rates(k,53) )

         mat(k,1438) = -( rxt(k,39) + rxt(k,40) + het_rates(k,54) )

         mat(k,453) = -( rxt(k,90) + het_rates(k,55) )

         mat(k,2105) = -( het_rates(k,56) )
         mat(k,243) = rxt(k,74)
         mat(k,107) = 4.000_r8*rxt(k,78)
         mat(k,151) = rxt(k,79)
         mat(k,122) = 2.000_r8*rxt(k,81)
         mat(k,161) = 2.000_r8*rxt(k,82)
         mat(k,126) = 2.000_r8*rxt(k,83)
         mat(k,166) = rxt(k,84)
         mat(k,130) = 2.000_r8*rxt(k,85)
         mat(k,133) = 3.000_r8*rxt(k,88)
         mat(k,376) = rxt(k,89)
         mat(k,184) = 2.000_r8*rxt(k,91)
         mat(k,101) = 2.000_r8*rxt(k,92)
         mat(k,2066) = rxt(k,93)
         mat(k,938) = rxt(k,94)
         mat(k,235) = rxt(k,97)
         mat(k,231) = rxt(k,100)
         mat(k,261) = rxt(k,101)
         mat(k,314) = rxt(k,102)
         mat(k,2190) = rxt(k,103)
         mat(k,811) = rxt(k,106)

         mat(k,183) = -( rxt(k,91) + het_rates(k,57) )

         mat(k,99) = -( rxt(k,92) + rxt(k,210) + het_rates(k,58) )

         mat(k,2065) = -( rxt(k,93) + het_rates(k,59) )
         mat(k,937) = rxt(k,95)
         mat(k,338) = rxt(k,107)
         mat(k,100) = 2.000_r8*rxt(k,210)

         mat(k,933) = -( rxt(k,94) + rxt(k,95) + rxt(k,533) + rxt(k,538) + rxt(k,544) &
                 + het_rates(k,60) )

         mat(k,5) = -( het_rates(k,61) )

         mat(k,1088) = -( het_rates(k,62) )
         mat(k,176) = 1.500_r8*rxt(k,22)
         mat(k,417) = .450_r8*rxt(k,23)
         mat(k,589) = .600_r8*rxt(k,26)
         mat(k,351) = rxt(k,27)
         mat(k,2025) = rxt(k,32) + rxt(k,33)
         mat(k,1043) = rxt(k,34)
         mat(k,1225) = rxt(k,36)
         mat(k,1436) = .380_r8*rxt(k,39)
         mat(k,835) = rxt(k,41)
         mat(k,1054) = rxt(k,43)
         mat(k,959) = 2.000_r8*rxt(k,44)
         mat(k,1096) = .330_r8*rxt(k,45)
         mat(k,1213) = 1.340_r8*rxt(k,51)
         mat(k,1336) = .700_r8*rxt(k,55)
         mat(k,206) = 1.500_r8*rxt(k,64)
         mat(k,653) = .250_r8*rxt(k,65)
         mat(k,999) = rxt(k,68)
         mat(k,1078) = 1.700_r8*rxt(k,69)
         mat(k,366) = rxt(k,110)

         mat(k,834) = -( rxt(k,41) + het_rates(k,63) )
         mat(k,617) = rxt(k,37)
         mat(k,1434) = .440_r8*rxt(k,39)
         mat(k,531) = .400_r8*rxt(k,60)
         mat(k,649) = rxt(k,65)
         mat(k,1075) = .800_r8*rxt(k,69)

         mat(k,253) = -( rxt(k,96) + het_rates(k,64) )
         mat(k,148) = rxt(k,79)
         mat(k,153) = rxt(k,80)
         mat(k,159) = rxt(k,82)
         mat(k,124) = 2.000_r8*rxt(k,83)
         mat(k,163) = 2.000_r8*rxt(k,84)
         mat(k,128) = rxt(k,85)
         mat(k,112) = 2.000_r8*rxt(k,98)
         mat(k,256) = rxt(k,101)
         mat(k,309) = rxt(k,102)

         mat(k,232) = -( rxt(k,97) + het_rates(k,65) )
         mat(k,120) = rxt(k,81)
         mat(k,158) = rxt(k,82)
         mat(k,228) = rxt(k,100)

         mat(k,200) = -( het_rates(k,66) )

         mat(k,303) = -( het_rates(k,67) )

         mat(k,6) = -( het_rates(k,68) )

         mat(k,7) = -( het_rates(k,69) )

         mat(k,8) = -( het_rates(k,70) )

         mat(k,9) = -( het_rates(k,71) )

         mat(k,10) = -( het_rates(k,72) )

         mat(k,11) = -( het_rates(k,73) )

         mat(k,12) = -( het_rates(k,74) )

         mat(k,13) = -( het_rates(k,75) )

         mat(k,14) = -( het_rates(k,76) )

         mat(k,15) = -( rxt(k,124) + het_rates(k,77) )

         mat(k,135) = -( rxt(k,42) + het_rates(k,78) )

         mat(k,879) = -( het_rates(k,79) )
         mat(k,154) = rxt(k,80)
         mat(k,164) = rxt(k,84)
         mat(k,254) = 2.000_r8*rxt(k,96)
         mat(k,233) = rxt(k,97)
         mat(k,292) = rxt(k,104)

      end do

      end subroutine linmat01

      subroutine linmat02( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
!       ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer,  intent(in)    ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(in)    ::  het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)
!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k


      do k = 1,avec_len

         mat(k,1053) = -( rxt(k,43) + het_rates(k,80) )
         mat(k,1095) = .330_r8*rxt(k,45)
         mat(k,651) = .250_r8*rxt(k,65)
         mat(k,276) = rxt(k,300)

         mat(k,958) = -( rxt(k,44) + rxt(k,480) + het_rates(k,81) )
         mat(k,400) = rxt(k,21)
         mat(k,416) = .130_r8*rxt(k,23)
         mat(k,342) = .700_r8*rxt(k,61)
         mat(k,666) = .600_r8*rxt(k,70)
         mat(k,712) = .340_r8*rxt(k,72)
         mat(k,429) = .170_r8*rxt(k,73)

         mat(k,1453) = -( rxt(k,138) + het_rates(k,82) )
         mat(k,2280) = rxt(k,2) + 2.000_r8*rxt(k,3)
         mat(k,2029) = 2.000_r8*rxt(k,33)
         mat(k,448) = rxt(k,38)
         mat(k,1439) = .330_r8*rxt(k,39) + rxt(k,40)
         mat(k,817) = rxt(k,99)
         mat(k,2179) = rxt(k,103)
         mat(k,293) = rxt(k,104)

         mat(k,1409) = -( het_rates(k,83) )
         mat(k,2277) = rxt(k,1)
         mat(k,2026) = rxt(k,32)
         mat(k,1437) = 1.440_r8*rxt(k,39)

         mat(k,111) = -( rxt(k,98) + het_rates(k,84) )

         mat(k,600) = -( rxt(k,4) + rxt(k,498) + het_rates(k,85) )
         mat(k,1632) = .050_r8*rxt(k,496)

         mat(k,138) = -( rxt(k,109) + het_rates(k,86) )

         mat(k,815) = -( rxt(k,99) + het_rates(k,87) )

         mat(k,227) = -( rxt(k,100) + het_rates(k,88) )

         mat(k,257) = -( rxt(k,101) + het_rates(k,89) )

         mat(k,310) = -( rxt(k,102) + het_rates(k,90) )

         mat(k,2192) = -( rxt(k,103) + het_rates(k,91) )

         mat(k,185) = -( het_rates(k,92) )

         mat(k,944) = -( het_rates(k,93) )
         mat(k,2024) = rxt(k,500)

         mat(k,291) = -( rxt(k,104) + het_rates(k,94) )

         mat(k,1471) = -( rxt(k,9) + het_rates(k,95) )
         mat(k,1101) = rxt(k,482)
         mat(k,582) = rxt(k,483)
         mat(k,551) = rxt(k,484)
         mat(k,280) = 2.000_r8*rxt(k,485) + 2.000_r8*rxt(k,494) + 2.000_r8*rxt(k,529) &
                      + 2.000_r8*rxt(k,532) + 2.000_r8*rxt(k,543)
         mat(k,382) = rxt(k,486)
         mat(k,1034) = rxt(k,487)
         mat(k,2224) = .500_r8*rxt(k,489)
         mat(k,1756) = rxt(k,490) + rxt(k,495)
         mat(k,388) = rxt(k,491)
         mat(k,248) = rxt(k,492)
         mat(k,625) = rxt(k,493)
         mat(k,565) = rxt(k,531) + rxt(k,534) + rxt(k,539)
         mat(k,934) = rxt(k,533) + rxt(k,538) + rxt(k,544)

         mat(k,421) = -( rxt(k,10) + rxt(k,11) + rxt(k,173) + het_rates(k,96) )

         mat(k,784) = -( rxt(k,105) + het_rates(k,97) )
         mat(k,563) = rxt(k,531) + rxt(k,534) + rxt(k,539)

         mat(k,808) = -( rxt(k,106) + het_rates(k,98) )
         mat(k,932) = rxt(k,533) + rxt(k,538) + rxt(k,544)

         mat(k,1097) = -( rxt(k,45) + rxt(k,482) + het_rates(k,99) )

         mat(k,236) = -( rxt(k,46) + het_rates(k,100) )
         mat(k,1278) = rxt(k,373)

         mat(k,1110) = -( rxt(k,47) + het_rates(k,101) )
         mat(k,1098) = .170_r8*rxt(k,45)

         mat(k,332) = -( het_rates(k,102) )

         mat(k,141) = -( het_rates(k,103) )

         mat(k,853) = -( het_rates(k,104) )

         mat(k,578) = -( rxt(k,483) + het_rates(k,105) )

         mat(k,546) = -( rxt(k,484) + het_rates(k,106) )

         mat(k,1198) = -( het_rates(k,107) )

         mat(k,391) = -( rxt(k,48) + het_rates(k,108) )

         mat(k,696) = -( rxt(k,49) + het_rates(k,109) )
         mat(k,392) = rxt(k,48)

         mat(k,80) = -( het_rates(k,110) )

         mat(k,1214) = -( rxt(k,50) + rxt(k,51) + het_rates(k,111) )
         mat(k,698) = .300_r8*rxt(k,49)

         mat(k,316) = -( het_rates(k,112) )

         mat(k,512) = -( rxt(k,52) + het_rates(k,113) )
         mat(k,637) = .800_r8*rxt(k,19)
         mat(k,671) = .800_r8*rxt(k,20)

         mat(k,321) = -( rxt(k,53) + het_rates(k,114) )

         mat(k,607) = -( rxt(k,54) + rxt(k,355) + het_rates(k,115) )

         mat(k,897) = -( het_rates(k,116) )

         mat(k,1340) = -( rxt(k,55) + het_rates(k,117) )
         mat(k,699) = .700_r8*rxt(k,49)

         mat(k,481) = -( rxt(k,155) + het_rates(k,118) )
         mat(k,1795) = rxt(k,15)

         mat(k,189) = -( rxt(k,12) + het_rates(k,119) )

         mat(k,279) = -( rxt(k,13) + rxt(k,14) + rxt(k,174) + rxt(k,485) + rxt(k,494) &
                      + rxt(k,529) + rxt(k,532) + rxt(k,543) + het_rates(k,120) )

         mat(k,379) = -( rxt(k,486) + het_rates(k,121) )

         mat(k,1030) = -( rxt(k,56) + rxt(k,487) + het_rates(k,122) )

         mat(k,16) = -( het_rates(k,123) )

         mat(k,17) = -( het_rates(k,124) )

         mat(k,18) = -( het_rates(k,125) )

         mat(k,102) = -( het_rates(k,126) )

         mat(k,19) = -( rxt(k,488) + het_rates(k,127) )

         mat(k,20) = -( rxt(k,547) + het_rates(k,128) )

         mat(k,21) = -( rxt(k,546) + het_rates(k,129) )

         mat(k,1853) = -( rxt(k,15) + het_rates(k,130) )
         mat(k,282) = rxt(k,14)
         mat(k,2229) = rxt(k,16) + .500_r8*rxt(k,489)
         mat(k,1761) = rxt(k,17)
         mat(k,485) = rxt(k,155)

         mat(k,2237) = -( rxt(k,16) + rxt(k,489) + het_rates(k,131) )
         mat(k,1481) = rxt(k,9)
         mat(k,425) = rxt(k,11) + rxt(k,173)
         mat(k,283) = rxt(k,13) + rxt(k,174)
         mat(k,1769) = rxt(k,18)
         mat(k,647) = rxt(k,19)
         mat(k,1107) = rxt(k,45)
         mat(k,396) = rxt(k,48)
         mat(k,615) = rxt(k,54) + rxt(k,355)
         mat(k,1040) = rxt(k,56)
         mat(k,842) = rxt(k,57)
         mat(k,390) = rxt(k,58)
         mat(k,249) = rxt(k,59)
         mat(k,537) = .600_r8*rxt(k,60) + rxt(k,306)
         mat(k,628) = rxt(k,66)
         mat(k,568) = rxt(k,76)
         mat(k,940) = rxt(k,95)
         mat(k,146) = rxt(k,430)

         mat(k,1760) = -( rxt(k,17) + rxt(k,18) + rxt(k,490) + rxt(k,495) &
                 + het_rates(k,132) )
         mat(k,423) = rxt(k,10)
         mat(k,281) = rxt(k,13) + rxt(k,14) + rxt(k,174)
         mat(k,534) = .400_r8*rxt(k,60)
         mat(k,567) = rxt(k,77)
         mat(k,936) = rxt(k,94)

         mat(k,837) = -( rxt(k,57) + het_rates(k,133) )

         mat(k,385) = -( rxt(k,58) + rxt(k,491) + het_rates(k,134) )

         mat(k,22) = -( het_rates(k,135) )

         mat(k,23) = -( het_rates(k,136) )

         mat(k,24) = -( het_rates(k,137) )

         mat(k,25) = -( het_rates(k,138) )

         mat(k,2269) = -( rxt(k,132) + het_rates(k,139) )
         mat(k,2295) = rxt(k,3)
         mat(k,2171) = rxt(k,8)
         mat(k,284) = rxt(k,14)
         mat(k,1862) = rxt(k,15)
         mat(k,2238) = rxt(k,16)
         mat(k,1770) = rxt(k,18)
         mat(k,1449) = .180_r8*rxt(k,39)
         mat(k,836) = rxt(k,41)
         mat(k,1506) = rxt(k,75)
         mat(k,2070) = rxt(k,93)
         mat(k,339) = rxt(k,107)
         mat(k,1248) = rxt(k,111) + rxt(k,473)
         mat(k,848) = rxt(k,112)
         mat(k,269) = rxt(k,113)
         mat(k,1547) = rxt(k,127) + rxt(k,128)
         mat(k,487) = rxt(k,155)
         mat(k,522) = rxt(k,466)

         mat(k,2168) = -( rxt(k,7) + rxt(k,8) + rxt(k,499) + het_rates(k,140) )
         mat(k,2266) = rxt(k,132)

         mat(k,26) = -( het_rates(k,141) )

         mat(k,335) = -( rxt(k,107) + het_rates(k,142) )

         mat(k,363) = -( rxt(k,110) + het_rates(k,143) )

         mat(k,247) = -( rxt(k,59) + rxt(k,492) + het_rates(k,144) )

         mat(k,530) = -( rxt(k,60) + rxt(k,306) + het_rates(k,145) )

         mat(k,144) = -( rxt(k,430) + het_rates(k,146) )

         mat(k,469) = -( het_rates(k,147) )
         mat(k,272) = rxt(k,31)

         mat(k,178) = -( het_rates(k,148) )

         mat(k,340) = -( rxt(k,61) + het_rates(k,149) )

         mat(k,27) = -( het_rates(k,150) )

         mat(k,28) = -( het_rates(k,151) )

         mat(k,538) = -( rxt(k,62) + het_rates(k,152) )

         mat(k,403) = -( rxt(k,63) + het_rates(k,153) )

         mat(k,517) = -( rxt(k,466) + het_rates(k,154) )
         mat(k,364) = rxt(k,110)
         mat(k,1235) = rxt(k,111)

         mat(k,29) = -( rxt(k,108) + het_rates(k,155) )

         mat(k,1237) = -( rxt(k,111) + rxt(k,473) + het_rates(k,156) )
         mat(k,845) = rxt(k,112)
         mat(k,518) = rxt(k,466)

         mat(k,844) = -( rxt(k,112) + het_rates(k,157) )
         mat(k,268) = rxt(k,113)
         mat(k,1236) = rxt(k,473)

         mat(k,267) = -( rxt(k,113) + het_rates(k,158) )
         mat(k,139) = rxt(k,109)

         mat(k,30) = -( het_rates(k,159) )

         mat(k,31) = -( het_rates(k,160) )

         mat(k,32) = -( het_rates(k,161) )

         mat(k,33) = -( rxt(k,114) + het_rates(k,162) )

         mat(k,34) = -( rxt(k,115) + het_rates(k,163) )

         mat(k,35) = -( rxt(k,116) + het_rates(k,164) )

         mat(k,36) = -( rxt(k,117) + het_rates(k,165) )

         mat(k,37) = -( rxt(k,118) + het_rates(k,166) )

         mat(k,38) = -( rxt(k,119) + het_rates(k,167) )

         mat(k,39) = -( rxt(k,120) + het_rates(k,168) )

         mat(k,40) = -( rxt(k,121) + het_rates(k,169) )

         mat(k,41) = -( rxt(k,122) + het_rates(k,170) )

         mat(k,42) = -( rxt(k,123) + het_rates(k,171) )

         mat(k,43) = -( het_rates(k,172) )
         mat(k,956) = rxt(k,480)

         mat(k,44) = -( het_rates(k,173) )

         mat(k,45) = -( het_rates(k,174) )

         mat(k,46) = -( het_rates(k,175) )

         mat(k,47) = -( het_rates(k,176) )

         mat(k,48) = -( rxt(k,548) + het_rates(k,177) )

         mat(k,54) = -( het_rates(k,178) )

         mat(k,205) = -( rxt(k,64) + het_rates(k,179) )

         mat(k,648) = -( rxt(k,65) + het_rates(k,180) )

      end do

      end subroutine linmat02

      subroutine linmat03( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
!       ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer,  intent(in)    ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(in)    ::  het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)
!----------------------------------------------
!       ... local variables
!----------------------------------------------
      integer :: k


      do k = 1,avec_len

         mat(k,623) = -( rxt(k,66) + rxt(k,493) + het_rates(k,181) )

         mat(k,488) = -( rxt(k,67) + het_rates(k,182) )

         mat(k,996) = -( rxt(k,68) + het_rates(k,183) )
         mat(k,386) = rxt(k,58)
         mat(k,624) = rxt(k,66)
         mat(k,490) = rxt(k,67)

         mat(k,1077) = -( rxt(k,69) + het_rates(k,184) )
         mat(k,652) = rxt(k,65)
         mat(k,998) = rxt(k,68)

         mat(k,661) = -( rxt(k,70) + het_rates(k,185) )

         mat(k,193) = -( het_rates(k,186) )

         mat(k,209) = -( rxt(k,71) + het_rates(k,187) )

         mat(k,214) = -( het_rates(k,188) )

         mat(k,707) = -( rxt(k,72) + het_rates(k,189) )

         mat(k,222) = -( het_rates(k,190) )

         mat(k,427) = -( rxt(k,73) + het_rates(k,191) )

         mat(k,524) = -( het_rates(k,194) )
         mat(k,145) = rxt(k,430)

         mat(k,921) = -( het_rates(k,195) )

         mat(k,60) = -( het_rates(k,196) )

         mat(k,497) = -( het_rates(k,197) )

         mat(k,66) = -( het_rates(k,198) )

         mat(k,435) = -( het_rates(k,199) )

         mat(k,826) = -( het_rates(k,200) )
         mat(k,514) = rxt(k,52)

         mat(k,793) = -( het_rates(k,201) )

         mat(k,631) = -( het_rates(k,202) )

         mat(k,1394) = -( het_rates(k,203) )
         mat(k,419) = .130_r8*rxt(k,23)
         mat(k,353) = rxt(k,27)
         mat(k,952) = rxt(k,35)
         mat(k,1227) = rxt(k,36)
         mat(k,1100) = .330_r8*rxt(k,45)
         mat(k,1112) = rxt(k,47)
         mat(k,1218) = 1.340_r8*rxt(k,50)
         mat(k,515) = rxt(k,52)
         mat(k,324) = rxt(k,53)
         mat(k,1342) = .300_r8*rxt(k,55)
         mat(k,839) = rxt(k,57)
         mat(k,532) = .600_r8*rxt(k,60) + rxt(k,306)
         mat(k,405) = rxt(k,63)
         mat(k,207) = .500_r8*rxt(k,64)
         mat(k,1080) = .650_r8*rxt(k,69)

         mat(k,1905) = -( het_rates(k,204) )
         mat(k,1048) = rxt(k,34)
         mat(k,954) = rxt(k,35)
         mat(k,620) = rxt(k,37)
         mat(k,1444) = rxt(k,40)
         mat(k,1348) = .300_r8*rxt(k,55)
         mat(k,535) = .400_r8*rxt(k,60)
         mat(k,595) = rxt(k,87)
         mat(k,374) = rxt(k,89)

         mat(k,765) = -( het_rates(k,205) )
         mat(k,301) = .600_r8*rxt(k,25)

         mat(k,570) = -( het_rates(k,206) )

         mat(k,275) = -( rxt(k,299) + rxt(k,300) + het_rates(k,207) )
         mat(k,136) = rxt(k,42)

         mat(k,720) = -( het_rates(k,208) )

         mat(k,2013) = -( rxt(k,481) + rxt(k,497) + het_rates(k,209) )
         mat(k,424) = rxt(k,11) + rxt(k,173)
         mat(k,645) = rxt(k,19)
         mat(k,680) = .900_r8*rxt(k,20)
         mat(k,402) = rxt(k,21)
         mat(k,177) = 1.500_r8*rxt(k,22)
         mat(k,420) = .560_r8*rxt(k,23)
         mat(k,468) = rxt(k,24)
         mat(k,302) = .600_r8*rxt(k,25)
         mat(k,590) = .600_r8*rxt(k,26)
         mat(k,354) = rxt(k,27)
         mat(k,331) = rxt(k,28)
         mat(k,349) = rxt(k,29)
         mat(k,413) = rxt(k,30)
         mat(k,1049) = rxt(k,34)
         mat(k,1231) = rxt(k,36)
         mat(k,1056) = 2.000_r8*rxt(k,43)
         mat(k,961) = 2.000_r8*rxt(k,44)
         mat(k,1105) = .670_r8*rxt(k,45)
         mat(k,240) = rxt(k,46)
         mat(k,1114) = rxt(k,47)
         mat(k,395) = rxt(k,48)
         mat(k,701) = rxt(k,49)
         mat(k,1220) = 1.340_r8*rxt(k,50) + .660_r8*rxt(k,51)
         mat(k,1038) = rxt(k,56)
         mat(k,344) = rxt(k,61)
         mat(k,543) = rxt(k,62)
         mat(k,208) = rxt(k,64)
         mat(k,656) = rxt(k,65)
         mat(k,627) = rxt(k,66)
         mat(k,493) = rxt(k,67)
         mat(k,1003) = rxt(k,68)
         mat(k,1083) = 1.200_r8*rxt(k,69)
         mat(k,670) = rxt(k,70)
         mat(k,717) = rxt(k,72)
         mat(k,432) = rxt(k,73)
         mat(k,1459) = rxt(k,138)
         mat(k,444) = rxt(k,269)
         mat(k,277) = rxt(k,299) + rxt(k,300)
         mat(k,1305) = rxt(k,373)

         mat(k,441) = -( rxt(k,269) + het_rates(k,210) )

         mat(k,1262) = -( het_rates(k,211) )

         mat(k,1295) = -( rxt(k,373) + het_rates(k,212) )

         mat(k,72) = -( het_rates(k,213) )

         mat(k,78) = -( het_rates(k,214) )

         mat(k,1319) = -( het_rates(k,215) )

         mat(k,727) = -( het_rates(k,216) )
         mat(k,467) = .600_r8*rxt(k,24)

         mat(k,1362) = -( het_rates(k,217) )
         mat(k,1217) = .660_r8*rxt(k,50)
         mat(k,610) = rxt(k,54) + rxt(k,355)

         mat(k,869) = -( het_rates(k,218) )
         mat(k,588) = .600_r8*rxt(k,26)

         mat(k,683) = -( het_rates(k,219) )

         mat(k,86) = -( het_rates(k,220) )

         mat(k,1063) = -( het_rates(k,221) )

         mat(k,1535) = -( rxt(k,127) + rxt(k,128) + het_rates(k,222) )
         mat(k,2283) = rxt(k,1)
         mat(k,2159) = rxt(k,7)
         mat(k,190) = rxt(k,12)

         mat(k,1700) = -( rxt(k,496) + het_rates(k,223) )
         mat(k,2284) = rxt(k,2)
         mat(k,601) = 2.000_r8*rxt(k,4)
         mat(k,1473) = rxt(k,9)
         mat(k,422) = rxt(k,10)
         mat(k,679) = rxt(k,20)
         mat(k,401) = rxt(k,21)
         mat(k,330) = rxt(k,28)
         mat(k,348) = rxt(k,29)
         mat(k,412) = rxt(k,30)
         mat(k,274) = rxt(k,31)
         mat(k,619) = rxt(k,37)
         mat(k,449) = rxt(k,38)
         mat(k,1442) = .330_r8*rxt(k,39)
         mat(k,137) = rxt(k,42)
         mat(k,239) = rxt(k,46)
         mat(k,700) = rxt(k,49)
         mat(k,325) = rxt(k,53)
         mat(k,389) = rxt(k,58)
         mat(k,343) = rxt(k,61)
         mat(k,542) = rxt(k,62)
         mat(k,406) = rxt(k,63)
         mat(k,655) = rxt(k,65)
         mat(k,492) = rxt(k,67)
         mat(k,669) = rxt(k,70)
         mat(k,211) = rxt(k,71)
         mat(k,716) = rxt(k,72)
         mat(k,431) = rxt(k,73)
         mat(k,787) = rxt(k,105)
         mat(k,809) = rxt(k,106)
         mat(k,2227) = .500_r8*rxt(k,489)

         mat(k,460) = -( het_rates(k,224) )

         mat(k,774) = -( het_rates(k,225) )

         mat(k,1161) = -( het_rates(k,226) )
         mat(k,1079) = .150_r8*rxt(k,69)

         mat(k,1124) = -( het_rates(k,227) )

         mat(k,1145) = -( het_rates(k,228) )

         mat(k,738) = -( het_rates(k,229) )

         mat(k,92) = -( het_rates(k,230) )

         mat(k,1178) = -( het_rates(k,231) )

         mat(k,754) = -( het_rates(k,232) )

         mat(k,98) = -( het_rates(k,233) )

         mat(k,505) = -( het_rates(k,234) )

         mat(k,2296) = -( rxt(k,1) + rxt(k,2) + rxt(k,3) + het_rates(k,235) )
         mat(k,1450) = .050_r8*rxt(k,39)
         mat(k,140) = rxt(k,109)
         mat(k,2021) = rxt(k,481) + .500_r8*rxt(k,497)
         mat(k,1712) = .450_r8*rxt(k,496)
         mat(k,606) = rxt(k,498)

      end do

      end subroutine linmat03

      subroutine linmat( avec_len, mat, y, rxt, het_rates )
!----------------------------------------------
!       ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
!       ... dummy arguments
!----------------------------------------------
      integer,  intent(in)    ::  avec_len
      real(r8), intent(in)    ::  y(veclen,gas_pcnst)
      real(r8), intent(in)    ::  rxt(veclen,rxntot)
      real(r8), intent(in)    ::  het_rates(veclen,gas_pcnst)
      real(r8), intent(inout) ::  mat(veclen,nzcnt)

      call linmat01( avec_len, mat, y, rxt, het_rates )
      call linmat02( avec_len, mat, y, rxt, het_rates )
      call linmat03( avec_len, mat, y, rxt, het_rates )

      end subroutine linmat

      end module mo_lin_matrix

