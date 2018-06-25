! Deallocate all the arrays inside body ibd
subroutine sm_deallocateBody ( ibd )
    use SolidMechanics_data, only :  sm_BodyInfo, sm_structure
    implicit none
    
    ! IO variables
    integer, intent(in) :: ibd
    
    ! Define internal variables
    type(sm_structure), pointer :: body
    integer :: irest

    body => sm_BodyInfo(ibd)
    
    if( allocated( body%x             ) ) deallocate( body%x             )
    if( allocated( body%y             ) ) deallocate( body%y             )
    if( allocated( body%z             ) ) deallocate( body%z             )
    if( allocated( body%IEN           ) ) deallocate( body%IEN           )
    if( allocated( body%ID            ) ) deallocate( body%ID            )
    if( allocated( body%LM            ) ) deallocate( body%LM            )
    if( allocated( body%MatDensity    ) ) deallocate( body%MatDensity    )
    if( allocated( body%YoungsModulus ) ) deallocate( body%YoungsModulus )
    if( allocated( body%PoissonsRatio ) ) deallocate( body%PoissonsRatio )
    if( allocated( body%MatType       ) ) deallocate( body%MatType       )
    if( allocated( body%eltype        ) ) deallocate( body%eltype        )
    if( allocated( body%qi            ) ) deallocate( body%qi            )
    if( allocated( body%K             ) ) deallocate( body%K             )
    if( allocated( body%M             ) ) deallocate( body%M             )
    if( allocated( body%Damp          ) ) deallocate( body%Damp          )
    if( allocated( body%qq_IA         ) ) deallocate( body%qq_IA         )
    if( allocated( body%qq_JA         ) ) deallocate( body%qq_JA         )
    if( allocated( body%Kqv           ) ) deallocate( body%Kqv           )
    if( allocated( body%Mqv           ) ) deallocate( body%Mqv           )
    if( allocated( body%Dampqv        ) ) deallocate( body%Dampqv        )
    if( allocated( body%qv_IA         ) ) deallocate( body%qv_IA         )
    if( allocated( body%qv_JA         ) ) deallocate( body%qv_JA         )
    if( allocated( body%fix_list_A    ) ) deallocate( body%fix_list_A    )
    if( allocated( body%LM_cs         ) ) deallocate( body%LM_cs         )
    if( allocated( body%ws_LM_cs      ) ) deallocate( body%ws_LM_cs      )
    if( allocated( body%Qs            ) ) deallocate( body%Qs            )
    if( allocated( body%ws_IEN        ) ) deallocate( body%ws_IEN        )
    if( allocated( body%ws_eltype     ) ) deallocate( body%ws_eltype     )
    if( allocated( body%Qsn           ) ) deallocate( body%Qsn           )
    if( allocated( body%Hs            ) ) deallocate( body%Hs            )
    if( allocated( body%Hsn           ) ) deallocate( body%Hsn           )
    if( allocated( body%Hs_pres       ) ) deallocate( body%Hs_pres       )
    if( allocated( body%Hs_visc       ) ) deallocate( body%Hs_visc       )
    if( allocated( body%Hsi_pres      ) ) deallocate( body%Hsi_pres      )
    if( allocated( body%Hsi_visc      ) ) deallocate( body%Hsi_visc      )
    if( allocated( body%dyn_rhs       ) ) deallocate( body%dyn_rhs       )
    if( allocated( body%qn            ) ) deallocate( body%qn            )
    if( allocated( body%qdn           ) ) deallocate( body%qdn           )
    if( allocated( body%qddn          ) ) deallocate( body%qddn          )
    if( allocated( body%qdi           ) ) deallocate( body%qdi           )
    if( allocated( body%qddi          ) ) deallocate( body%qddi          )
    if( allocated( body%qms           ) ) deallocate( body%qms           )
    if( allocated( body%qdms          ) ) deallocate( body%qdms          )
    if( allocated( body%qddms         ) ) deallocate( body%qddms         )
    if( allocated( body%ey            ) ) deallocate( body%ey            )
    if( allocated( body%edy           ) ) deallocate( body%edy           )
    if( allocated( body%restraints    ) ) deallocate( body%restraints    )
    if( allocated( body%restraints_surf ) ) then
      do irest = 1,SIZE(body%restraints_surf)
         if( allocated( body%restraints_surf(irest)%node_list ) ) then
            deallocate( body%restraints_surf(irest)%node_list )
         end if
      end do
    end if
    
end
