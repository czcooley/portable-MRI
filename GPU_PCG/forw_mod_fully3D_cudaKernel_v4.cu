


// this version of the kernel supports PE
// this is the pixel dependent field variation from rot to rot expressed in radians


__device__ void multiply_2complex( float a_re,float a_im,
				   float b_re,float b_im,
				   float *c_re,float *c_im )
{
  *c_re = a_re*b_re - a_im*b_im;
  *c_im = a_re*b_im + a_im*b_re;

}




__global__ void forw_mod_fully3D_cudaKernel_v1( float *y_re,float *y_im,
						float *x_re,float *x_im,
						float *times,float *fieldmaps,
						float *b1m_re,float *b1m_im,
						float *b1p_re,float *b1p_im,
						float *b1p_echoe_exp,
						float *freqs,float *freq_weights,
                        float *PE_field,
                        int nfreqs,
						int nsamples,int ncoils,int nrots,
						int nechoes,int npix )
{
  
  int index = threadIdx.x + blockIdx.x * blockDim.x;

  if( index < nsamples*ncoils*nrots*nechoes*nfreqs ){
    
    int tmp1 = index % (nsamples*ncoils*nrots*nechoes);
    int tmp2 = tmp1 % (nsamples*ncoils*nrots);
    int tmp3 = tmp2 % (nsamples*ncoils);
    int tmp4 = tmp3 % nsamples;

    int sample = tmp4;
    int coil   = (tmp3 - tmp4 ) / nsamples;
    int rot    = (tmp2 - tmp3 ) / (nsamples*ncoils);
    int echoe  = (tmp1 - tmp2 ) / (nsamples*ncoils*nrots);
    int freq   = (index - tmp1) / (nsamples*ncoils*nrots*nechoes);
    
    y_re[index] = 0.0f;
    y_im[index] = 0.0f;

    for( int i=0;i<npix;i++ ){
      
      // EXP
      float angle = -times[sample] * ( fieldmaps[i + rot*npix] + freqs[freq]  )   +    PE_field[ i + rot*npix ];
      float cos_angle = cosf(angle);
      float sin_angle = sinf(angle);
      
      // B1-
      float b1m_re_ = b1m_re[ i + coil*npix + rot*npix*ncoils ];
      float b1m_im_ = b1m_im[ i + coil*npix + rot*npix*ncoils ];
      
      // B1+
      float b1p_re_ = b1p_re[ i + rot*npix ];
      float b1p_im_ = b1p_im[ i + rot*npix ];

      angle = atan2f( b1p_im_,b1p_re_ ) * b1p_echoe_exp[echoe];
      float mag = sqrtf( b1p_re_*b1p_re_ +  b1p_im_*b1p_im_ );

      b1p_re_ = mag * cosf(angle);
      b1p_im_ = mag * sinf(angle);
      

      // compute matrix element for this 
      float A,B,C,D,E,F;
      multiply_2complex( cos_angle,sin_angle,
			 b1m_re_,b1m_im_,
			 &A,&B);

      multiply_2complex( A,B,
			 b1p_re_,b1p_im_,
			 &C,&D);      

      // matrix multiplication
      multiply_2complex( C,D, 
			 x_re[i],x_im[i],
			 &E,&F );

      y_re[index] += E * freq_weights[freq];
      y_im[index] += F * freq_weights[freq];

      
    }
  

  }

}

































