#ifdef __cplusplus
extern "C" {
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#ifdef __cplusplus
}
#endif

MODULE = Math::WalshTransform	PACKAGE = Math::WalshTransform
PROTOTYPES: ENABLE

void
xs_fht (n, ...) 
CODE:
{
	register int i;
	double *mr;
	int j, k, l, nk, nl;

	/* Or, declare n as an arg */
	unsigned int n;
	n = (int) SvIV(ST(0));
	/* perlapi SvNV for a double, SvIV for an int, SvPV for a char* etc */
	
	/* New() is perlapi's equivalent to malloc(). See perlclib.pod */
	New(0, mr, n, double);
	for (i=0; i<n; i++) mr[i] = (double) SvNV(ST(i+1));

   k = 1;
   l = n;
   while (1) {
      i = -1; l = l/2;
      for (nl=1; nl<=l; nl++) {
         for (nk=1; nk<=k; nk++) {
            i++; j = i+k;
				/* fprintf (stderr, "xs_fht: nk=%d  nl=%d i=%d j=%d mr[i]=%g\n",
				 * nk, nl, i, j, mr[i]); */
            mr[i] = (mr[i] + mr[j])/2;
            mr[j] =  mr[i] - mr[j];
         }
         i = j;
      }
      k = 2*k;
      if (k >= n) { break; }
   }
   if (k == n) {
		for (i=0; i<n; i++) ST(i) = sv_2mortal(newSVnv(mr[i]));
		Safefree(mr);
		XSRETURN(n);
   } else {
		fprintf (stderr, "fht: n should be a power of 2, but was %d\n", n);
		XSRETURN_EMPTY;
   }
}

void
xs_fhtinv (n, ...) 
CODE:
{
	register int i;
	double *mr;
	int j, k, l, nk, nl;

	/* Or, declare n as an arg */
	unsigned int n;
	n = (int) SvIV(ST(0));
	/* perlapi SvNV for a double, SvIV for an int, SvPV for a char* etc */
	
	/* New() is perlapi's equivalent to malloc(). See perlclib.pod */
	New(0, mr, n, double);
	for (i=0; i<n; i++) mr[i] = (double) SvNV(ST(i+1));

   k = 1;
   l = n;
   while (1) {
      i = -1; l = l/2;
      for (nl=1; nl<=l; nl++) {
         for (nk=1; nk<=k; nk++) {
            i++; j = i+k;
				/* fprintf (stderr, "xs_fht: nk=%d  nl=%d i=%d j=%d mr[i]=%g\n",
				 * nk, nl, i, j, mr[i]); */
            mr[i] = mr[i] + mr[j];
            mr[j] = mr[i] - 2.0*mr[j];
         }
         i = j;
      }
      k = 2*k;
      if (k >= n) { break; }
   }
   if (k == n) {
		for (i=0; i<n; i++) ST(i) = sv_2mortal(newSVnv(mr[i]));
		Safefree(mr);
		XSRETURN(n);
   } else {
		fprintf (stderr, "fhtinv: n should be a power of 2, but was %d\n", n);
		XSRETURN_EMPTY;
   }
}

void
xs_fwt (n, ...) 
CODE:
{
	register int i;
	double *mr;
	double *nr;
	int k, l, m, nh, nk, nl, tmp, alternate, kp1, kh;

	unsigned int n;
	n = (int) SvIV(ST(0)); /* perlapi SvNV=double, SvIV=int, SvPV=char* etc */
	
	New(0, mr, n, double); /* perlapi's malloc(). See perlclib.pod */
	for (i=0; i<n; i++) mr[i] = (double) SvNV(ST(i+1));
	New(0, nr, n, double);

	m = 0; tmp = 1;
	while (1) { if (tmp>=n) break; tmp<<=1; m++; }
	alternate = m & 1;

	if (alternate) {
		for (k=0; k<n; k+=2) {
			kp1 = k+1;
			mr[k]   = 0.5 * (mr[k] + mr[kp1]);
			mr[kp1] =  mr[k] - mr[kp1];
		}
	} else {
		for (k=0; k<n; k+=2) {
			kp1 = k+1;
			nr[k]   = 0.5 * (mr[k] + mr[kp1]);
			nr[kp1] =  nr[k] - mr[kp1];
		}
	}

	k = 1; nh = n/2;
	while (1) {
		kh = k; k <<= 1; kp1 = k+1; if (kp1>n) break;
		nh = nh/2; l = 0; i = 0; alternate = !alternate;
		for (nl=1; nl<=nh; nl++) {
			for (nk=1; nk<=kh; nk++) {
				if (alternate) {
					mr[l]   = 0.5 * (nr[i] + nr[i+k]);
					mr[l+1] = mr[l] - nr[i+k];
					mr[l+2] = 0.5 * (nr[i+1] - nr[i+kp1]);
					mr[l+3] = mr[l+2] + nr[i+kp1];
				} else {
					nr[l]   = 0.5 * (mr[i] + mr[i+k]);
					nr[l+1] = nr[l] - mr[i+k];
					nr[l+2] = 0.5 * (mr[i+1] - mr[i+kp1]);
					nr[l+3] = nr[l+2] + mr[i+kp1];
				}
				l = l+4; i = i+2;
			}
			i = i+k;
		}
	}
	Safefree(nr);
	for (i=0; i<n; i++) ST(i) = sv_2mortal(newSVnv(mr[i]));
	Safefree(mr);
	XSRETURN(n);
}

void
xs_fwtinv (n, ...) 
CODE:
{
	register int i;
	double *mr;
	double *nr;
	int k, l, m, nh, nk, nl, tmp, alternate, kp1, kh;

	unsigned int n;
	n = (int) SvIV(ST(0)); /* perlapi SvNV=double, SvIV=int, SvPV=char* etc */
	
	New(0, mr, n, double); /* perlapi's malloc(). See perlclib.pod */
	for (i=0; i<n; i++) mr[i] = (double) SvNV(ST(i+1));
	New(0, nr, n, double);

	m = 0; tmp = 1;
	while (1) { if (tmp>=n) break; tmp<<=1; m++; }
	alternate = m & 1;

	if (alternate) {
		for (k=0; k<n; k+=2) {
			kp1 = k+1;
			mr[k]   = mr[k] + mr[kp1];
			mr[kp1] = mr[k] - mr[kp1] - mr[kp1];
		}
	} else {
		for (k=0; k<n; k+=2) {
			kp1 = k+1;
			nr[k]   = mr[k] + mr[kp1];
			nr[kp1] = mr[k] - mr[kp1];
		}
	}

	k = 1; nh = n/2;
	while (1) {
		kh = k; k <<= 1; kp1 = k+1; if (kp1>n) break;
		nh = nh/2; l = 0; i = 0; alternate = !alternate;
		for (nl=1; nl<=nh; nl++) {
			for (nk=1; nk<=kh; nk++) {
				if (alternate) {
					mr[l]   = nr[i]   + nr[i+k];
					mr[l+1] = nr[i]   - nr[i+k];
					mr[l+2] = nr[i+1] - nr[i+kp1];
					mr[l+3] = nr[i+1] + nr[i+kp1];
				} else {
					nr[l]   = mr[i]   + mr[i+k];
					nr[l+1] = mr[i]   - mr[i+k];
					nr[l+2] = mr[i+1] - mr[i+kp1];
					nr[l+3] = mr[i+1] + mr[i+kp1];
				}
				l = l+4; i = i+2;
			}
			i = i+k;
		}
	}
	Safefree(nr);
	for (i=0; i<n; i++) ST(i) = sv_2mortal(newSVnv(mr[i]));
	Safefree(mr);
	XSRETURN(n);
}
