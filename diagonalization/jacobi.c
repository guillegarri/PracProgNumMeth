#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>

int rotations;

void reset_rotations() {
  rotations = 0;
}

int get_rotations(){
  return rotations;
}

int jacobi(gsl_matrix * A, gsl_vector * e, gsl_matrix * V){
  int changed , sweeps=0, n=A->size1;


  for (int i = 0; i < n; i++) {
    gsl_vector_set(e,i,gsl_matrix_get(A,i,i));
  }
  gsl_matrix_set_identity(V);


  do {changed=0;
    sweeps++;
    int p,q;

    for (p = 0; p < n; p++) {
      for (q = p+1; q < n; q++) {
        double App = gsl_vector_get(e,p);
        double Aqq = gsl_vector_get(e,q);
        double Apq = gsl_matrix_get(A,p,q);

        double phi = atan2(2*Apq, Aqq-App)/2.0;
        rotations++;

        double cosphi = cos(phi);
        double sinphi = sin(phi);

        double AppNew = cosphi*cosphi*App - 2*sinphi*cosphi*Apq + sinphi*sinphi*Aqq;
        double AqqNew = sinphi*sinphi*App + 2*sinphi*cosphi*Apq + cosphi*cosphi*Aqq;

        if (AppNew != App || AqqNew != Aqq) {
          changed = 1;
          gsl_vector_set(e,p,AppNew);
          gsl_vector_set(e,q,AqqNew);
          gsl_matrix_set(A,p,q,0.0);

          for (int i = 0; i < p; i++) {
            double Aip = gsl_matrix_get(A,i,p);
            double Aiq = gsl_matrix_get(A,i,q);

            gsl_matrix_set(A,i,p, cosphi*Aip - sinphi*Aiq);
            gsl_matrix_set(A,i,q, cosphi*Aiq + sinphi*Aip);
          }
          for (int i = p+1; i < q; i++) {
            double Api = gsl_matrix_get(A,p,i);
            double Aiq = gsl_matrix_get(A,i,q);

            gsl_matrix_set(A,p,i, cosphi*Api - sinphi*Aiq);
            gsl_matrix_set(A,i,q, cosphi*Aiq + sinphi*Api);
          }
          for (int i = q+1; i < n; i++) {
            double Api = gsl_matrix_get(A,p,i);
            double Aqi = gsl_matrix_get(A,q,i);

            gsl_matrix_set(A,p,i, cosphi*Api - sinphi*Aqi);
            gsl_matrix_set(A,q,i, cosphi*Aqi + sinphi*Api);
          }
          for (int i = 0; i < n; i++) {
            double Vip = gsl_matrix_get(V,i,p);
            double Viq = gsl_matrix_get(V,i,q);

            gsl_matrix_set(V,i,p, cosphi*Vip - sinphi*Viq);
            gsl_matrix_set(V,i,q, cosphi*Viq + sinphi*Vip);
          }
        }
      }
    }
  } while(changed != 0);
  return sweeps;
}

int jacobiN_LtH(gsl_matrix * A, gsl_vector * e, gsl_matrix * V, int numvals) {
  int changed , sweeps=0, n=A->size1;


  for (int i = 0; i < n; i++) {
    gsl_vector_set(e,i,gsl_matrix_get(A,i,i));
  }
  gsl_matrix_set_identity(V);

    int p,q;

    for (p = 0; p < numvals; p++) {
      do {changed=0;
        sweeps++;
      for (q = p+1; q < n; q++) {
        double App = gsl_vector_get(e,p);
        double Aqq = gsl_vector_get(e,q);
        double Apq = gsl_matrix_get(A,p,q);

        double phi = atan2(2*Apq, Aqq-App)/2.0;
        rotations++;

        double cosphi = cos(phi);
        double sinphi = sin(phi);

        double AppNew = cosphi*cosphi*App - 2*sinphi*cosphi*Apq + sinphi*sinphi*Aqq;
        double AqqNew = sinphi*sinphi*App + 2*sinphi*cosphi*Apq + cosphi*cosphi*Aqq;

        if (AppNew != App || AqqNew != Aqq) {
          changed = 1;
          gsl_vector_set(e,p,AppNew);
          gsl_vector_set(e,q,AqqNew);
          gsl_matrix_set(A,p,q,0.0);

          for (int i = 0; i < p; i++) {
            double Aip = gsl_matrix_get(A,i,p);
            double Aiq = gsl_matrix_get(A,i,q);

            gsl_matrix_set(A,i,p, cosphi*Aip - sinphi*Aiq);
            gsl_matrix_set(A,i,q, cosphi*Aiq + sinphi*Aip);
          }
          for (int i = p+1; i < q; i++) {
            double Api = gsl_matrix_get(A,p,i);
            double Aiq = gsl_matrix_get(A,i,q);

            gsl_matrix_set(A,p,i, cosphi*Api - sinphi*Aiq);
            gsl_matrix_set(A,i,q, cosphi*Aiq + sinphi*Api);
          }
          for (int i = q+1; i < n; i++) {
            double Api = gsl_matrix_get(A,p,i);
            double Aqi = gsl_matrix_get(A,q,i);

            gsl_matrix_set(A,p,i, cosphi*Api - sinphi*Aqi);
            gsl_matrix_set(A,q,i, cosphi*Aqi + sinphi*Api);
          }
          for (int i = 0; i < n; i++) {
            double Vip = gsl_matrix_get(V,i,p);
            double Viq = gsl_matrix_get(V,i,q);

            gsl_matrix_set(V,i,p, cosphi*Vip - sinphi*Viq);
            gsl_matrix_set(V,i,q, cosphi*Viq + sinphi*Vip);
          }
        }
      }
    }
    while(changed != 0);
  }

  for (int i = numvals; i < n; i++) {
    gsl_vector_set(e,i,0);
  }
  return sweeps;
}

int jacobiN_HtL(gsl_matrix * A, gsl_vector * e, gsl_matrix * V, int numvals) {
  int changed , sweeps=0, n=A->size1;


  for (int i = 0; i < n; i++) {
    gsl_vector_set(e,i,gsl_matrix_get(A,i,i));
  }
  gsl_matrix_set_identity(V);

    int p,q;

    for (p = 0; p < numvals; p++) {
      do {changed=0;
        sweeps++;
      for (q = p+1; q < n; q++) {
        double App = gsl_vector_get(e,p);
        double Aqq = gsl_vector_get(e,q);
        double Apq = gsl_matrix_get(A,p,q);

        double phi = atan2(2*Apq, Aqq-App)/2.0 + M_PI/2;
        rotations++;

        double cosphi = cos(phi);
        double sinphi = sin(phi);

        double AppNew = cosphi*cosphi*App - 2*sinphi*cosphi*Apq + sinphi*sinphi*Aqq;
        double AqqNew = sinphi*sinphi*App + 2*sinphi*cosphi*Apq + cosphi*cosphi*Aqq;

        if (AppNew != App || AqqNew != Aqq) {
          changed = 1;
          gsl_vector_set(e,p,AppNew);
          gsl_vector_set(e,q,AqqNew);
          gsl_matrix_set(A,p,q,0.0);

          for (int i = 0; i < p; i++) {
            double Aip = gsl_matrix_get(A,i,p);
            double Aiq = gsl_matrix_get(A,i,q);

            gsl_matrix_set(A,i,p, cosphi*Aip - sinphi*Aiq);
            gsl_matrix_set(A,i,q, cosphi*Aiq + sinphi*Aip);
          }
          for (int i = p+1; i < q; i++) {
            double Api = gsl_matrix_get(A,p,i);
            double Aiq = gsl_matrix_get(A,i,q);

            gsl_matrix_set(A,p,i, cosphi*Api - sinphi*Aiq);
            gsl_matrix_set(A,i,q, cosphi*Aiq + sinphi*Api);
          }
          for (int i = q+1; i < n; i++) {
            double Api = gsl_matrix_get(A,p,i);
            double Aqi = gsl_matrix_get(A,q,i);

            gsl_matrix_set(A,p,i, cosphi*Api - sinphi*Aqi);
            gsl_matrix_set(A,q,i, cosphi*Aqi + sinphi*Api);
          }
          for (int i = 0; i < n; i++) {
            double Vip = gsl_matrix_get(V,i,p);
            double Viq = gsl_matrix_get(V,i,q);

            gsl_matrix_set(V,i,p, cosphi*Vip - sinphi*Viq);
            gsl_matrix_set(V,i,q, cosphi*Viq + sinphi*Vip);
          }
        }
      }
    }
    while(changed != 0);
  }

  for (int i = numvals; i < n; i++) {
    gsl_vector_set(e,i,0);
  }
  return sweeps;
}
