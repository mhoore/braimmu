#include "math.h"
#include "cstring"
#include "nii2pov.h"

#define SECOND 1.0/2.0
#define THIRD 1.0/3.0

#define CAM_ANGLE 33.0
#define ZOOM_DIST 170.0
#define SCALE_BAR 50.0

#define TRANS 0.8
#define FILTER 0.0

#define LOGSCALE 0
#define SMALL 1.e-10

/* ----------------------------------------------------------------------*/
N2P::N2P(int narg, char **arg)
{
  if (narg < 7) {
    printf ("ERROR: arg1 = input, arg2 = output,\
            arg3 = which_data_to_show (integer),\
            arg4 = minval, arg5 = maxval, arg6 = section (integer) \n");
    exit (EXIT_FAILURE);
  }

  if (LOGSCALE)
    printf ("Warning: logscale colorbar is turned on. \n");

  me = atoi(arg[3]);
  minval = atof(arg[4]);
  maxval = atof(arg[5]);
  sec_id = atoi(arg[6]);

  read_image(arg[1]);
  mri_boundaries(nim);
  mri_topology(nim);

  sprintf(fname,arg[2]);
  fw = fopen(fname,"w");

  settings();
  file_write();

  fclose(fw);

}

/* ----------------------------------------------------------------------*/
N2P::~N2P()
{
  destroy(data);
  destroy(x);
}

/* ----------------------------------------------------------------------*/
void N2P::read_image(char *arg) {
  nim = nifti_image_read(arg,1);

  if(!nim) {
    printf("Error: cannot open nifti image %s \n", arg);
    exit(1);
  }

  printf("##################### \n");
  printf("NIFTI image %s is read. \n", arg);
  printf("NIFTI image properties: ");
  printf("ndim = %i \n", nim->ndim);
  for (int i=1; i<8; i++)
    printf("dim[%i] = %i, pixdim[i] = %g \n",
           i,nim->dim[i],i,nim->pixdim[i]);
  printf("nvox = %lli \n", nim->nvox);
  printf("nbyper = %i \n", nim->nbyper);
  printf("datatype = %i \n", nim->datatype);

  printf("calmin = %g, calmax = %g \n", nim->cal_min, nim->cal_max);
  printf("toffset = %g \n", nim->toffset);

  printf("xyz_units = %i, time_units = %i \n", nim->xyz_units, nim->time_units);
  printf("nifti_type = %i \n", nim->nifti_type);

  printf("intent_code = %i \n", nim->intent_code);

  printf("description: %s \n", nim->descrip);
  printf("##################### \n");

}

/* ----------------------------------------------------------------------*/
int N2P::mri_boundaries(nifti_image *nim) {
  int i;

  double conver_fac = 1.0;
  if (nim->xyz_units == NIFTI_UNITS_METER)
    conver_fac = 1.e6;
  else if (nim->xyz_units == NIFTI_UNITS_MM)
    conver_fac = 1.e3;
  else if (nim->xyz_units == NIFTI_UNITS_MICRON)
    conver_fac = 1.0;

  for (i=0; i<3; i++) {
    vlen[i] = nim->pixdim[i+1] * conver_fac;
    lbox[i] = nim->dim[i+1] * nim->pixdim[i+1];
    lbox[i] *= conver_fac;
    boxlo[i] = 0.0;
    boxhi[i] = lbox[i];
    nv[i] = nim->dim[i+1];
  }

  return 1;

}

/* ----------------------------------------------------------------------*/
void N2P::settings () {
  int i;
  ifstream fr;
  string line;

  double light[2][3], dum;

  fr.open("settings.pov");
  fprintf(fw, "/// POV-RAY image \n");
  do {
    getline(fr, line);
    fprintf(fw, "%s\n", line.c_str());
  } while(!fr.eof());
  fr.close();

  //initial camera positioning
  cam[0] = 1;
  cam[1] = 0;
  cam[2] = 0;

  up[0] = 0;
  up[1] = 0;
  up[2] = 1;

  right[0] = 0;
  right[1] = 1;
  right[2] = 0;

  look[0] = 0.5 * lbox[0];
  look[1] = 0.5 * lbox[1];
  look[2] = 0.5 * lbox[2];

  // up normalization
  dum = 1.0/sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);
  for(i = 0; i < 3; i++)
    up[i] *= dum;

  for (i=0; i<3; i++) {
    light[0][i] = (cam[i] + up[i] + right[i]) * 2.0 * lbox[i] + look[i];
    light[1][i] = (cam[i] + up[i] - right[i]) * 2.0 * lbox[i] + look[i];
  }

  //camera position and lighting
  fprintf(fw, "//camera and lights \n");
  fprintf(fw, "camera{ orthographic angle  %g \n", CAM_ANGLE);
  fprintf(fw, "  location <350000, 200000, 300000> \n");
//  fprintf(fw, "  location <%g, %g, %g> \n",
//          0.5 * lbox[0] + cam[0] * 2.0 * lbox[0],
//          0.5 * lbox[1] + cam[1] * 2.0 * lbox[1],
//          0.5 * lbox[2] + cam[2] * 2.0 * lbox[2]);

  fprintf(fw, "  look_at  <%g, %g, %g> \n",
          look[0],
          look[1],
          look[2]);

  fprintf(fw, "  right <%g,%g,%g>*image_width/image_height \n", right[0], right[1], right[2]);
  fprintf(fw, "  up <%g,%g,%g> }\n", up[0], up[1], up[2]);

  fprintf(fw, "light_source { <%g,%g,%g> color White} \n", light[0][0], light[0][1], light[0][2]);
  fprintf(fw, "light_source { <%g,%g,%g> color White} \n", light[1][0], light[1][1], light[1][2]);
  fprintf(fw, "\n\n");

}

/* ----------------------------------------------------------------------*/
void N2P::file_write () {
  tagint c;
  int i,j,k;
  double x0,y0,z0,x1,y1,z1;
  double scbar0[3],scbar1[3], down, lbar;

  double grayscale, *rgb, coef;

  create(rgb,3,"rgb");

  // if not specified, find min and max value from the data
  if (minval == 0 && maxval == 0) {
    minval = data[0][me];
    maxval = data[0][me];

    for (c=0; c<nvoxel; c++) {
      if (minval > data[c][me])
        minval = data[c][me];
      if (maxval < data[c][me])
        maxval = data[c][me];
    }
  }

  if (LOGSCALE) {
    maxval = log10(maxval + SMALL);
    minval = log10(minval + SMALL);
  }

  coef = 2.0 / (maxval - minval);

  if (nim->dim[0] == 3) {
    c = 0;
    for (k=0; k<nim->dim[3]; k++) {
      for (j=0; j<nim->dim[2]; j++) {
        for (i=0; i<nim->dim[1]; i++) {

          if (i != sec_id) {
            c++;
            continue;
          }

          x0 = x[c][0] - 0.5*vlen[0];
          x1 = x[c][0] + 0.5*vlen[0];

          y0 = x[c][1] - 0.5*vlen[1];
          y1 = x[c][1] + 0.5*vlen[1];

          z0 = x[c][2] - 0.5*vlen[2];
          z1 = x[c][2] + 0.5*vlen[2];

          if (LOGSCALE)
            grayscale = (log10(data[c][0] + SMALL) - minval) * coef - 1.0;
          else
            grayscale = (data[c][0] - minval) * coef - 1.0;
          jet(grayscale, rgb);

          fprintf(fw, "box{<%g, %g, %g>, <%g, %g, %g> \
                  texture{pigment{color rgbft <%g, %g, %g, %g, %g>} } }\n",
                  x0, y0, z0, x1, y1, z1, rgb[0], rgb[1], rgb[2], TRANS, FILTER);

          c++;
        }
      }
    }
  }

  else if (nim->dim[0] == 5) {

    // if not specified, find min and max value from the data
    if (minval == 0 && maxval == 0) {
      minval = data[0][me];
      maxval = data[0][me];

      for (c=0; c<nvoxel; c++) {
        if (minval > data[c][me])
          minval = data[c][me];
        if (maxval < data[c][me])
          maxval = data[c][me];
      }
    }

    coef = 2.0 / (maxval - minval);

    c = 0;
    for (k=0; k<nim->dim[3]; k++) {
      for (j=0; j<nim->dim[2]; j++) {
        for (i=0; i<nim->dim[1]; i++) {

          //if (i != sec_id || data[c][dsize-1] <= 0) {
          if (data[c][dsize-1] <= 0) {
            c++;
            continue;
          }

          if ((x[c][0] - 0.4 * lbox[0]) >= 0 &&
              (x[c][1] - 0.4 * lbox[1]) >= 0 &&
              (x[c][2] - 0.4 * lbox[2]) >= 0 ) {
            c++;
            continue;
          }

          if ((x[c][2] - 0.65 * lbox[2]) >= 0 ) {
            c++;
            continue;
          }

          x0 = x[c][0] - 0.5*vlen[0];
          x1 = x[c][0] + 0.5*vlen[0];

          y0 = x[c][1] - 0.5*vlen[1];
          y1 = x[c][1] + 0.5*vlen[1];

          z0 = x[c][2] - 0.5*vlen[2];
          z1 = x[c][2] + 0.5*vlen[2];

          grayscale = (data[c][me] - minval) * coef - 1.0;
          jet(grayscale, rgb);

/*
          if (data[c][me] == 30 || data[c][me] == 17) {
            rgb[0] = 0; rgb[1] = 0; rgb[2] = 0;
          }
          else if (data[c][me] == 210 || data[c][me] == 211) {
            rgb[0] = 0; rgb[1] = 0; rgb[2] = 255.0/256;
          }
          else if (data[c][me] == 83 || data[c][me] == 59) {
            rgb[0] = 255.0/256; rgb[1] = 0; rgb[2] = 0;
          }
          else if (data[c][me] == 218 || data[c][me] == 219) {
            rgb[0] = 0; rgb[1] = 255.0/256; rgb[2] = 0;
          }
          else if (data[c][me] == 57 || data[c][me] == 105) {
            rgb[0] = 0; rgb[1] = 255.0/256; rgb[2] = 255.0/256;
          }
          else if (data[c][me] == 6 || data[c][me] == 2) {
            rgb[0] = 255.0/256; rgb[1] = 0; rgb[2] = 255.0/256;
          }
          else if (data[c][me] == 73 || data[c][me] == 45) {
            rgb[0] = 128.0/256; rgb[1] = 128.0/256; rgb[2] = 0;
          }
          else if (data[c][me] == 8 || data[c][me] == 4) {
            rgb[0] = 144.0/256; rgb[1] = 238.0/256; rgb[2] = 144.0/256;
          }
          else if (data[c][me] == 67 || data[c][me] == 76) {
            rgb[0] = 178.0/256; rgb[1] = 34.0/256; rgb[2] = 34.0/256;
          }
          else if (data[c][me] == 39 || data[c][me] == 53) {
            rgb[0] = 255.0/256; rgb[1] = 128.0/256; rgb[2] = 0;
          }
          else if (data[c][me] == 14 || data[c][me] == 16) {
            rgb[0] = 70.0/256; rgb[1] = 130.0/256; rgb[2] = 180.0/256;
          }
          else if (data[c][me] == 102 || data[c][me] == 203) {
            rgb[0] = 255.0/256; rgb[1] = 105.0/256; rgb[2] = 180.0/256;
          }
          else if (data[c][me] == 33 || data[c][me] == 23) {
            rgb[0] = 112.0/256; rgb[1] = 128.0/256; rgb[2] = 144.0/256;
          }
          else if (data[c][me] == 12 || data[c][me] == 11) {
            rgb[0] = 138.0/256; rgb[1] = 43.0/256; rgb[2] = 226.0/256;
          }
          else {
            c++;
            continue;
          }
*/
          fprintf(fw, "box{<%g, %g, %g>, <%g, %g, %g> \
                  texture{pigment{color rgbft <%g, %g, %g, %g, %g>} } }\n",
                  x0, y0, z0, x1, y1, z1, rgb[0], rgb[1], rgb[2], TRANS, FILTER);

          c++;
        }
      }
    }
  }

  //scale bar
  down = -ZOOM_DIST * tan((0.5*CAM_ANGLE)/180*M_PI) + 2.0;
  lbar = 0.5 * SCALE_BAR;
  for (i=0; i<3; i++ ) {
    scbar0[i] = down * up[i] + right[i] * lbar + look[i];
    scbar1[i] = down * up[i] - right[i] * lbar + look[i];
  }

  //fprintf(f_write, "cylinder { <%g,%g,%g>, <%g,%g,%g>, 0.5 texture{pigment{color rgbft <0,0,0,0,0>}}}\n",scbar0[0],scbar0[1],scbar0[2],scbar1[0],scbar1[1],scbar1[2]);

  destroy(rgb);

}

// ----------------------------------------------------------------------
void N2P::mri_topology(nifti_image *nim) {
  int i,j,k,h;
  tagint c;

  if (nim->datatype == DT_UINT8)
    ptr8 = (uint8_t *) nim->data;
  else if (nim->datatype == DT_INT16)
    ptr16 = (int16_t *) nim->data;
  else if (nim->datatype == DT_INT32)
    ptr32 = (int32_t *) nim->data;
  else if (nim->datatype == DT_FLOAT32)
    ptrf = (float *) nim->data;
  else if (nim->datatype == DT_FLOAT64)
    ptrd = (double *) nim->data;
  else {
    printf("Error: nifti file data type cannot be read. datatype=%i . \n", nim->datatype);
    exit(1);
  }

  if (nim->dim[0] == 3)
    dsize = 1;
  else if (nim->dim[0] == 5)
    dsize = nim->dim[5];
  else {
    printf("Error: cannot read nifti file data. nim->dim[0]= %i. \n", nim->dim[0]);
    exit(1);
  }

  nvoxel = nim->dim[3] * nim->dim[2] * nim->dim[1];

  create(data,nvoxel,dsize,"data");
  create(x,nvoxel,3,"x");

  double dumx,dumy,dumz;

  double conver_fac = 1.0;
  if (nim->xyz_units == NIFTI_UNITS_METER)
    conver_fac = 1.e6;
  else if (nim->xyz_units == NIFTI_UNITS_MM)
    conver_fac = 1.e3;
  else if (nim->xyz_units == NIFTI_UNITS_MICRON)
    conver_fac = 1.0;

  // go through the nifti_image data and find the corresponding voxel
  if (nim->ndim == 3) {
    c = 0;
    for (k=0; k<nim->dim[3]; k++) {
      dumz = nim->pixdim[3] * k * conver_fac;

      for (j=0; j<nim->dim[2]; j++) {
        dumy = nim->pixdim[2] * j * conver_fac;

        for (i=0; i<nim->dim[1]; i++) {
          dumx = nim->pixdim[1] * i * conver_fac;

          if (nim->datatype == DT_UINT8)
            data[c][0] = (double) ptr8[c];
          else if (nim->datatype == DT_INT16)
            data[c][0] = (double) ptr16[c];
          else if (nim->datatype == DT_INT32)
            data[c][0] = (double) ptr32[c];
          else if (nim->datatype == DT_FLOAT32)
            data[c][0] = (double) ptrf[c];
          else if (nim->datatype == DT_FLOAT64)
            data[c][0] = (double) ptrd[c];

          x[c][0] = dumx;
          x[c][1] = dumy;
          x[c][2] = dumz;

          c++;
        }
      }
    }
  }

  else if (nim->ndim == 5) {
    tagint c0 = 0;
    for (h=0; h<nim->dim[5]; h++) {
      c = 0;
      for (k=0; k<nim->dim[3]; k++) {
        dumz = nim->pixdim[3] * k * conver_fac;

        for (j=0; j<nim->dim[2]; j++) {
          dumy = nim->pixdim[2] * j * conver_fac;

          for (i=0; i<nim->dim[1]; i++) {
            dumx = nim->pixdim[1] * i * conver_fac;

            // if it is in the partition
            if (nim->datatype == DT_UINT8)
              data[c][h] += (double) ptr8[c0];
            else if (nim->datatype == DT_INT16)
              data[c][h] += (double) ptr16[c0];
            else if (nim->datatype == DT_INT32)
              data[c][h] += (double) ptr32[c0];
            else if (nim->datatype == DT_FLOAT32)
              data[c][h] += (double) ptrf[c0];
            else if (nim->datatype == DT_FLOAT64)
              data[c][h] += (double) ptrd[c0];

            x[c][0] = dumx;
            x[c][1] = dumy;
            x[c][2] = dumz;

            c++;
            c0++;
          }
        }
      }
    }

    string str = nim->descrip;
    darg.clear();

    istringstream buf(str);

    for(string word; buf >> word;)
    darg.push_back(word);

    int dsize0 = darg.size();

    if (dsize0 != dsize){
      printf("Error: mri file contents do not match its description. \n");
      exit(1);
    }

  }

}

/* ----------------------------------------------------------------------*/
void N2P::sfree(void *ptr) {
  if (ptr == NULL) return;
  free(ptr);
}

/* ----------------------------------------------------------------------*/
void *N2P::smalloc(size_t nbytes, const char *name) {
  if (nbytes == 0) return NULL;

  void *ptr = malloc(nbytes);
  if (ptr == NULL) {
    printf("Failed to allocate %lld bytes for array %s. \n", nbytes,name);
    exit(1);
  }
  return ptr;
}

/* ----------------------------------------------------------------------*/
double N2P::interpolate(double val, double y0, double x0, double y1, double x1) {
  return (val-x0)*(y1-y0)/(x1-x0) + y0;
}

/* ----------------------------------------------------------------------*/
void N2P::jet(double grayscale, double *rgb) {
  // blue
  if ( grayscale < -0.33 ) rgb[2] = 1.0;
  else if ( grayscale < 0.33 ) rgb[2] = interpolate(grayscale, 1.0, -0.33, 0.0, 0.33 );
  else rgb[2] = 0.0;

  // green
  if ( grayscale < -1.0 ) rgb[1] = 0.0; // unexpected grayscale value
  if  ( grayscale < -0.33 ) rgb[1] = interpolate( grayscale, 0.0, -1.0, 1.0, -0.33 );
  else if ( grayscale < 0.33 ) rgb[1] = 1.0;
  else if ( grayscale <= 1.0 ) rgb[1] = interpolate( grayscale, 1.0, 0.33, 0.0, 1.0 );
  else rgb[1] = 1.0; // unexpected grayscale value

  // red
  if ( grayscale < -0.33 ) rgb[0] = 0.0;
  else if ( grayscale < 0.33 ) rgb[0] = interpolate( grayscale, 0.0, -0.33, 1.0, 0.33 );
  else rgb[0] = 1.0;

}

/* ----------------------------------------------------------------------*/
int main(int argc, char **argv) {
  N2P fdata(argc,argv);
  return 0;
}
