#include "fdebug.h"

module actuator_line_source

    use fldebug
    use spud
    use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi
    use actuator_line_model

    implicit none
    real,save :: meshFactor, dragFactor, chordFactor
    real, allocatable :: Sx(:),Sy(:),Sz(:),Sc(:),Se(:),Su(:),Sv(:),Sw(:),SFX(:),SFY(:),SFZ(:)
    integer :: NSource
 
    public get_locations, get_forces, set_vel, initialize_actuator_source
    
contains
    
    subroutine initialize_actuator_source

    implicit none
    integer :: counter,itur,iblade,ielem,ial
    
    ewrite(1,*) 'entering initialize_source_terms'  
    !> Get Source term parameters
    call get_option("/turbine_models/actuator_line_model/source_term_parameters/meshFactor",meshFactor,default=2.00)
    call get_option("/turbine_models/actuator_line_model/source_term_parameters/dragFactor",dragFactor,default=1.00)
    call get_option("/turbine_models/actuator_line_model/source_term_parameters/chordFactor",chordFactor,default=0.25)
    
    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            
            ! Blades
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                end do
            end do
            
            ! Tower 
            if(turbine(itur)%has_tower) then
                do ielem=1,Turbine(itur)%Tower%Nelem
                counter=counter+1
                end do
            endif
            
            ! Hub 
            if(turbine(itur)%has_hub) then
                do ielem=1,Turbine(itur)%hub%Nelem
                counter=counter+1
                end do
            endif

        end do
    endif
    
    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
            end do
        end do
    endif
    NSource=counter
    allocate(Sx(NSource),Sy(NSource),Sz(NSource),Sc(Nsource),Su(NSource),Sv(NSource),Sw(NSource),Se(NSource),Sfx(NSource),Sfy(NSource),Sfz(NSource))

    ewrite(1,*) 'exiting initialize_source_terms'

    end subroutine initialize_actuator_source

    subroutine get_locations
    
    implicit none
    integer :: counter,itur,iblade,ielem,ial

    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                Sx(counter)=Turbine(itur)%Blade(iblade)%PEX(ielem)
                Sy(counter)=Turbine(itur)%Blade(iblade)%PEY(ielem)
                Sz(counter)=Turbine(itur)%Blade(iblade)%PEZ(ielem)
                Sc(counter)=Turbine(itur)%Blade(iblade)%EC(ielem)
                end do
            end do
                !Tower 
                if(turbine(itur)%has_tower) then
                do ielem=1,Turbine(itur)%Tower%Nelem
                counter=counter+1
                Sx(counter)=Turbine(itur)%Tower%PEX(ielem)
                Sy(counter)=Turbine(itur)%Tower%PEY(ielem)
                Sz(counter)=Turbine(itur)%Tower%PEZ(ielem)
                Sc(counter)=Turbine(itur)%Tower%EC(ielem) 
                end do
                endif
                ! Hub 
                if(turbine(itur)%has_hub) then
                do ielem=1,Turbine(itur)%hub%Nelem
                counter=counter+1
                Sx(counter)=Turbine(itur)%hub%PEX(ielem)
                Sy(counter)=Turbine(itur)%hub%PEY(ielem)
                Sz(counter)=Turbine(itur)%hub%PEZ(ielem)
                Sc(counter)=Turbine(itur)%hub%EC(ielem) 
                end do
                endif
                
        end do
    endif
    
    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
                Sx(counter)=actuatorline(ial)%PEX(ielem)
                Sy(counter)=actuatorline(ial)%PEY(ielem)
                Sz(counter)=actuatorline(ial)%PEZ(ielem)
                Sc(counter)=actuatorline(ial)%EC(ielem)
            end do
        end do
    endif
  
    end subroutine get_locations
    
    subroutine set_vel
    
    implicit none
    integer :: counter,itur,iblade,ielem,ial
    
    ewrite(1,*) 'Entering set_vel'
    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            ! Blades
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                Turbine(itur)%Blade(iblade)%EVX(ielem)=Su(counter)
                Turbine(itur)%Blade(iblade)%EVY(ielem)=Sv(counter)
                Turbine(itur)%Blade(iblade)%EVZ(ielem)=Sw(counter)
                Turbine(itur)%Blade(iblade)%Eepsilon(ielem)=Se(counter) 
                end do
            end do
            ! Tower
            if(turbine(itur)%has_tower) then
                do ielem=1,Turbine(itur)%Tower%Nelem
                counter=counter+1
                Turbine(itur)%Tower%EVX(ielem)=Su(counter)
                Turbine(itur)%Tower%EVY(ielem)=Sv(counter)
                Turbine(itur)%Tower%EVZ(ielem)=Sw(counter)
                Turbine(itur)%Tower%Eepsilon(ielem)=Se(counter) 
                end do
            endif
            
            ! Hub
            if(turbine(itur)%has_hub) then
                do ielem=1,Turbine(itur)%hub%Nelem
                counter=counter+1
                Turbine(itur)%hub%EVX(ielem)=Su(counter)
                Turbine(itur)%hub%EVY(ielem)=Sv(counter)
                Turbine(itur)%hub%EVZ(ielem)=Sw(counter)
                Turbine(itur)%hub%Eepsilon(ielem)=Se(counter) 
                end do
            endif

        end do
    endif
    
    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
                actuatorline(ial)%EVX(ielem)=Su(counter)
                actuatorline(ial)%EVY(ielem)=Sv(counter)
                actuatorline(ial)%EVZ(ielem)=Sw(counter)
                actuatorline(ial)%Eepsilon(ielem)=Se(counter)
            end do
        end do
    endif

    ewrite(1,*) 'Exiting set_vel'
  
    end subroutine set_vel
    
    subroutine get_forces
    
    implicit none
    integer :: counter,itur,iblade,ielem,ial

    counter=0
    if (Ntur>0) then
        do itur=1,Ntur
            !Blade
            do iblade=1,Turbine(itur)%Nblades
                do ielem=1,Turbine(itur)%Blade(iblade)%Nelem
                counter=counter+1
                Sfx(counter)=Turbine(itur)%Blade(iblade)%EFX(ielem)
                Sfy(counter)=Turbine(itur)%Blade(iblade)%EFY(ielem)
                Sfz(counter)=Turbine(itur)%Blade(iblade)%EFZ(ielem)
                end do
            end do
            
            !Tower 
            if(turbine(itur)%has_tower) then
                do ielem=1,Turbine(itur)%Tower%Nelem
                counter=counter+1
                Sfx(counter)=Turbine(itur)%Tower%EFX(ielem)
                Sfy(counter)=Turbine(itur)%Tower%EFY(ielem)
                Sfz(counter)=Turbine(itur)%Tower%EFZ(ielem)
                end do
            endif
            
            ! Hub 
            if(turbine(itur)%has_hub) then
                do ielem=1,Turbine(itur)%hub%Nelem
                counter=counter+1
                Sfx(counter)=Turbine(itur)%hub%EFX(ielem)
                Sfy(counter)=Turbine(itur)%hub%EFY(ielem)
                Sfz(counter)=Turbine(itur)%hub%EFZ(ielem)
                end do
            endif

        end do
    endif
    
    if (Nal>0) then
        do ial=1,Nal
            do ielem=1,actuatorline(ial)%NElem
                counter=counter+1
                Sfx(counter)=actuatorline(ial)%EFX(ielem)
                Sfy(counter)=actuatorline(ial)%EFY(ielem)
                Sfz(counter)=actuatorline(ial)%EFZ(ielem)
            end do
        end do
    endif
  
    end subroutine get_forces

end module actuator_line_source
