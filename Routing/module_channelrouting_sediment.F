module module_channelrouting_sediment

implicit none

contains

!!DX--------------------------------------------------------------------------------------------
!!DX Begin Subroutine ChRt_sed()
!!DX--------------------------------------------------------------------------------------------

subroutine ChRt_sed(qlink_chan_sed,Bw_chan_sed,ChSSlp_chan_sed,hfrom_chan_sed,hto_chan_sed,&
                    zfrom_chan_sed,zto_chan_sed,Chlen_chan_sed,DTCT_chan_sed,&
                    jfrom_chan_sed,kfrom_chan_sed,jto_chan_sed,kto_chan_sed)

use module_namelist, only: nlst_rt

use module_RT_data, only: rt_domain

use module_overlandrouting_sediment


!!DX--------------------------------------------------------------------------------------------
!!DX Define input variables
!!DX--------------------------------------------------------------------------------------------

real,intent(in)::qlink_chan_sed,Bw_chan_sed,ChSSlp_chan_sed,hfrom_chan_sed,hto_chan_sed,&

                 zfrom_chan_sed,zto_chan_sed,Chlen_chan_sed,DTCT_chan_sed

integer,intent(in)::jfrom_chan_sed,kfrom_chan_sed,jto_chan_sed,kto_chan_sed




!!DX--------------------------------------------------------------------------------------------
!!DX Define local variables
!!DX--------------------------------------------------------------------------------------------


real ::v_chan_sed,supply_chan_sed,Rhfrom_chan_sed,sf_chan_sed,zto_chan_sed_use,&

       test_z,ADV_chan_sed,Trtot_chan

real :: Trtot_chan_sed(SizeFr),CS_chan_sed(SizeFr),SPG,gravity,dia(SizeFr),&

        EHcap(SizeFr),Suscap_chan_sed(SizeFr),ExcessCap_chan_sed(SizeFr),&

        Dep_chan_sed(SizeFr)

!! SPG is the specific gravity of the sediment, here we assume as 2.65, should
!depend upon the obs data, gravity is the gravity acceleartion,here is 9.81

real,dimension(ix_sed,jx_sed) :: SedFlux_chan

real,dimension(ix_sed,jx_sed,SizeFr) :: SedFlux_size_chan

integer :: size_chan_sed



!!DX--------------------------------------------------------------------------------------------
!!DX Initialize local variables
!!DX--------------------------------------------------------------------------------------------

supply_chan_sed=0.

SPG=2.65

gravity=9.81

dia(1)=0.00035 !! Sand unit:m

dia(2)=0.000016 !! Silt

dia(3)=0.000001 !! Clay

Trtot_chan=0.

SedFlux_chan=0.

SedFlux_size_chan=0.

! Test the bad elevation data

test_z=(zfrom_chan_sed-zto_chan_sed)/(Chlen_chan_sed)

if(test_z .lt. 0.002) then

zto_chan_sed_use=zfrom_chan_sed-Chlen_chan_sed*0.002
        
else
          
zto_chan_sed_use=zto_chan_sed
        
endif

sf_chan_sed=((zfrom_chan_sed-zto_chan_sed_use)+(hfrom_chan_sed-hto_chan_sed))/Chlen_chan_sed !! Caculate the friction slope


v_chan_sed=qlink_chan_sed/((Bw_chan_sed+hfrom_chan_sed/ChSSlp_chan_sed)*hfrom_chan_sed)

!!DX--------------------------------------------------------------------------------------------
!!DX Begin sediment transport 
!!DX--------------------------------------------------------------------------------------------


do size_chan_sed=1,SizeFr

supply_chan_sed=Sus(jfrom_chan_sed,kfrom_chan_sed,size_chan_sed)+Dep(jfrom_chan_sed,kfrom_chan_sed,size_chan_sed)

Trtot_chan_sed(size_chan_sed)=0.

if (supply_chan_sed .ne.0) then

!!DX--------------------------------------------------------------------------------------------
!!DX Calculate the sediment transport capacity 
!!DX--------------------------------------------------------------------------------------------


Rhfrom_chan_sed=((Bw_chan_sed+hfrom_chan_sed/ChSSlp_chan_sed)*hfrom_chan_sed)/&
                (Bw_chan_sed+2*hfrom_chan_sed*sqrt(1+1/(ChSSlp_chan_sed*ChSSlp_chan_sed))) !!Hydraulic Radius

CS_chan_sed(size_chan_sed)=0.05*(SPG/(SPG-1))*((v_chan_sed*sf_chan_sed)/(sqrt((SPG-1)*gravity*dia(size_chan_sed))))*&

                           (sqrt((Rhfrom_chan_sed*sf_chan_sed)/((SPG-1)*dia(size_chan_sed))))

!! dia(size_chan_sed) is the diameter of the corresponding size fraction

EHcap(size_chan_sed)=qlink_chan_sed*CS_chan_sed(size_chan_sed)*DTCT_chan_sed

ADV_chan_sed=min(v_chan_sed*DTCT_chan_sed/Chlen_chan_sed,1.0)

Suscap_chan_sed(size_chan_sed)=Sus(jfrom_chan_sed,kfrom_chan_sed,size_chan_sed)*ADV_chan_sed

CALL TRANSFERED(size_chan_sed,Sus,Suscap_chan_sed,Sus,jfrom_chan_sed,kfrom_chan_sed,jto_chan_sed,kto_chan_sed)

!! Compute the excess transport capacity

ExcessCap_chan_sed(size_chan_sed)=max(0.0,EHcap(size_chan_sed)-Suscap_chan_sed(size_chan_sed))

Dep_chan_sed(size_chan_sed)=Dep(jfrom_chan_sed,kfrom_chan_sed,size_chan_sed)*ADV_chan_sed

Dep_chan_sed(size_chan_sed)=min(ExcessCap_chan_sed(size_chan_sed),Dep_chan_sed(size_chan_sed))

CALL TRANSFERED(size_chan_sed,Dep,Dep_chan_sed,Sus,jfrom_chan_sed,kfrom_chan_sed,jto_chan_sed,kto_chan_sed)

Trtot_chan_sed(size_chan_sed)=Suscap_chan_sed(size_chan_sed)+Dep_chan_sed(size_chan_sed)

SedFlux_chan(jfrom_chan_sed,kfrom_chan_sed)=Trtot_chan_sed(size_chan_sed)+SedFlux_chan(jfrom_chan_sed,kfrom_chan_sed)

SedFlux_size_chan(jfrom_chan_sed,kfrom_chan_sed,size_chan_sed)=Trtot_chan_sed(size_chan_sed)/DTCT_chan_sed

end if

end do

end subroutine ChRt_sed
end module module_channelrouting_sediment
