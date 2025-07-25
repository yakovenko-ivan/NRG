module kind_parameters
   implicit none
   public

   !> Single precision real numbers, 6 digits, range 10?�7 to 10�7-1; 32 bits
   integer, parameter :: sp = selected_real_kind(6, 37)
   !> Double precision real numbers, 15 digits, range 10?��7 to 10��7-1; 64 bits
   integer, parameter :: dp = selected_real_kind(15, 307)
   !> Quadruple precision real numbers, 33 digits, range 10?4?�� to 104?��-1; 128 bits
   integer, parameter :: qp = selected_real_kind(33, 4931)

   !> Char length for integers, range -27 to 27-1; 8 bits
   integer, parameter :: i1 = selected_int_kind(2)
   !> Short length for integers, range -2�5 to 2�5-1; 16 bits
   integer, parameter :: i2 = selected_int_kind(4)
   !> Length of default integers, range -2�� to 2��-1; 32 bits
   integer, parameter :: i4 = selected_int_kind(9)
   !> Long length for integers, range -26� to 26�-1; 64 bits
   integer, parameter :: i8 = selected_int_kind(18)

end module kind_parameters
