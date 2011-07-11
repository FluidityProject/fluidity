  subroutine addref_REFCOUNT_TYPE(object)
#ifdef _OPENMP
    use omp_lib
#endif
    !!< Increment the reference count of object creating a new reference
    !!< counter if needed.
    type(REFCOUNT_TYPE), intent(inout), target :: object
    integer, save :: id = 0

    if (associated(object%refcount)) then
    !$OMP CRITICAL
       ! Reference count already exists, just increment it.
       object%refcount%count=object%refcount%count+1
    !$OMP END CRITICAL
       
    else
    !$OMP CRITICAL
       id = id + 1
       object%refcount=>new_refcount("REFCOUNT_TYPE", object%name)
       object%refcount%id = id
    !$OMP END CRITICAL
    end if
    
  end subroutine addref_REFCOUNT_TYPE
  
  subroutine incref_REFCOUNT_TYPE(object)
    !!< Increment the reference count of object. If there are no references
    !!< then error.
    type(REFCOUNT_TYPE), intent(in), target :: object
    integer, pointer :: ptr !! Dummy pointer to evade compilers which
    !! don't understand the rules for intent.

    if (.not.associated(object%refcount)) then
       FLAbort ("Attempt to incref REFCOUNT_TYPE "//trim(object%name)//" which has no references")
    end if
       
    ! Reference count already exists, just increment it.
    !$OMP CRITICAL
    ptr=>object%refcount%count
    ptr=ptr+1
    !$OMP END CRITICAL

  end subroutine incref_REFCOUNT_TYPE  
  
  subroutine decref_REFCOUNT_TYPE(object)
    !!< Decrement the reference count on object. If the reference count drops
    !!< to 0 deallocate the refcount as a hint to the calling routine that
    !!< the object can safely be deallocated.
    type(REFCOUNT_TYPE), intent(inout) :: object
    
    if (.not.associated(object%refcount)) then
       ! No refcount. Just exit
       return
    end if

    !$OMP CRITICAL
    object%refcount%count=object%refcount%count-1
    !$OMP END CRITICAL

    !$OMP CRITICAL
    if (object%refcount%count<=0) then

       if (object%refcount%count<0) then
          ! Warn for negative reference count
          ewrite(0,'(a, i0)') "Reference count of &
               &REFCOUNT_TYPE "//trim(object%name)//&
               " is ", object%refcount%count
          FLAbort("that should never happen.")
       end if

       object%refcount%prev%next=>object%refcount%next
       if (associated(object%refcount%next)) then
          object%refcount%next%prev=>object%refcount%prev
       end if
       
       deallocate(object%refcount)
       
    end if
    !$OMP END CRITICAL

  end subroutine decref_REFCOUNT_TYPE

  pure function has_references_REFCOUNT_TYPE(object) result (has_references)
    !!< Return true if there are any references to object
    type(REFCOUNT_TYPE), intent(in) :: object
    logical :: has_references
    
    has_references=associated(object%refcount)

  end function has_references_REFCOUNT_TYPE

