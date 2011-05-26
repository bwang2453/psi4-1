#include <stdio.h>
#include <string.h>
#include "libint.h"
#include "hrr_header.h"

extern void vrr_order_d0hd(Libint_t*, prim_data*);

  /* Computes quartets of (d0|hd) integrals */

REALTYPE *hrr_order_d0hd(Libint_t *Libint, int num_prim_comb)
{
 prim_data *Data = Libint->PrimQuartet;
 REALTYPE *int_stack = Libint->int_stack;
 int i;

 Libint->vrr_classes[2][5] = int_stack + 0;
 Libint->vrr_classes[2][6] = int_stack + 126;
 Libint->vrr_classes[2][7] = int_stack + 294;
 memset(int_stack,0,510*sizeof(REALTYPE));

 Libint->vrr_stack = int_stack + 510;
 for(i=0;i<num_prim_comb;i++) {
   vrr_order_d0hd(Libint, Data);
   Data++;
 }
 /*--- compute (d0|hp) ---*/
 hrr3_build_hp(Libint->CD,int_stack+510,int_stack+126,int_stack+0,6);
 /*--- compute (d0|ip) ---*/
 hrr3_build_ip(Libint->CD,int_stack+888,int_stack+294,int_stack+126,6);
 /*--- compute (d0|hd) ---*/
 hrr3_build_hd(Libint->CD,int_stack+1392,int_stack+888,int_stack+510,6);
 return int_stack+1392;}