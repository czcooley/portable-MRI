



// this version of the kernel supports PE
// this is the pixel dependent field variation from rot to rot expressed in radians



__device__ void multiply_2complex( float a_re,float a_im,
				   float b_re,float b_im,
				   float *c_re,float *c_im )
{
  *c_re = a_re*b_re - a_im*b_im;
  *c_im = a_re*b_im + a_im*b_re;

}




__global__ void forw_mod_fully3D_cudaKernel_TRANSP_v1( float *x_re,float *x_im,
						       float *y_re,float *y_im,
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

  if( index < npix ){
    
    x_re[index] = 0.0f;
    x_im[index] = 0.0f;

    for( int freq=0;freq<nfreqs;freq++ ){

      for( int sample=0;sample<nsamples;sample++ ){
	
	for( int rot=0;rot<nrots;rot++ ){
	  
	  // EXP
	  float angle = -times[sample] * ( fieldmaps[index + rot*npix] + freqs[freq] )   +    PE_field[ index + rot*npix ];;
	  float cos_angle = cosf(angle);
	  float sin_angle = sinf(angle);
	  
	  for( int coil=0;coil<ncoils;coil++ ){
	    
	    // B1-
	    float b1m_re_ = b1m_re[ index + coil*npix + rot*npix*ncoils ];
	    float b1m_im_ = b1m_im[ index + coil*npix + rot*npix*ncoils ];
	    
	    for( int echoe=0;echoe<nechoes;echoe++ ){      
	      
	      // B1+
	      float b1p_re_ = b1p_re[ index + rot*npix ];
	      float b1p_im_ = b1p_im[ index + rot*npix ];

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
	      int data_ind = sample + coil*nsamples + rot*nsamples*ncoils + echoe*nsamples*ncoils*nrots + freq*nsamples*ncoils*nrots*nechoes;
	      
	      multiply_2complex( C,-D, 
				 y_re[data_ind],y_im[data_ind],
				 &E,&F );
	      
	      x_re[index] += E * freq_weights[freq];
	      x_im[index] += F * freq_weights[freq];  

	    }

	  }

	}

      }
      
    }
    
    
  }

}

































