#include <math.h>

// Simple edge detection
SDoublePlane gradient(const SDoublePlane &in)
{
  SDoublePlane result(in.rows(), in.cols());

  for(int i=1; i<in.rows()-1; i++)
    for(int j=1; j<in.cols()-1; j++)
      {
	 double x = in[i-1][j-1] * -1 + in[i-1][j+1] * 1 + in[i][j-1] * -2 + in[i][j+1] * 2 + in[i+1][j-1] * -1 + in[i+1][j+1] * 1;
	 double y = in[i-1][j-1] * -1 + in[i-1][j] * -2 + in[i-1][j+1] * -1 + in[i+1][j-1] * 1 + in[i+1][j] * 2 + in[i+1][j+1] * 1;
	 result[i][j] = sqrt(x*x + y*y);
      }
    
  return result;
}

// This class takes care of the unary potentials (phi_i's) for you! 
// It learns parameters from the training data, and then can apply the
//  learned models to compute the unary potentials on new images.
//
class UnaryPotential
{
public:
  UnaryPotential(const vector<GTImage> &training_data, int _template_size=15) : template_size(_template_size)
  {
    // learn a little 31x31 template for each part. The template consists of a 
    //  estimate of mean gradient magnitue and variance in gradient magnitude 
    //  at each position.
    //

    // initialize models
    for(int i=0; i<6; i++)
      {
	means[i] = SDoublePlane(template_size*2+1, template_size*2+1);
	vars[i] = SDoublePlane(template_size*2+1, template_size*2+1);
      }  

    // train
    for(int img=0; img<training_data.size(); img++)
      {
	const GTImage &gt = training_data[img];

	// load the image
	cerr << "Loading training image " << gt.image_filename << "..." << endl;
	SDoublePlane image1 = SImageIO::read_png_file(gt.image_filename.c_str());

	// compute image gradients
	SDoublePlane grad = gradient(image1);

	// now update the model for this image, for each of 6 parts
	for(int ii=0; ii<6; ii++)
	  {
	    for(int i=-template_size; i<=template_size; i++)
	      for(int j=-template_size; j<=template_size; j++)
		{
		  double delta = grad[int(gt.coords[ii].row)+i][int(gt.coords[ii].col)+j] - means[ii][i+template_size][j+template_size];
		  means[ii][i+template_size][j+template_size] = means[ii][i+template_size][j+template_size] + delta/(img+1);
		  vars[ii][i+template_size][j+template_size] += delta * (grad[int(gt.coords[ii].row)+i][int(gt.coords[ii].col)+j] - means[ii][i+template_size][j+template_size]);
		}
	  }
      }

    for(int ii=0; ii<6; ii++)
      for(int i=-template_size; i<=template_size; i++)
	for(int j=-template_size; j<=template_size; j++)
	  vars[ii][i+template_size][j+template_size] /= (training_data.size());
  }

  // computes log potentials given a new image
  vector<SDoublePlane> compute_unary_logpotential(const string &image_filename) 
  {
    cerr << "Loading training image " << image_filename << "..." << endl;
    SDoublePlane image1 = SImageIO::read_png_file(image_filename.c_str());

    // compute image gradients
    SDoublePlane grad = gradient(image1);
    vector<SDoublePlane> unary_potentials(6);

    cerr << "Computing unary potentials..." << endl;
    for(int ii=0; ii<6; ii++)
      {
	SDoublePlane res(image1.rows(), image1.cols());
	for(int i=0; i<image1.rows(); i++)
	  for(int j=0; j<image1.cols(); j++)
	    res[i][j] = -1e100;

	for(int i=template_size; i<image1.rows()-template_size; i++)
	  for(int j=template_size; j<image1.cols()-template_size; j++)
	    {
	      res[i][j] = 0;
	      for(int i2=-template_size; i2<=template_size; i2++)
		for(int j2=-template_size; j2<=template_size; j2++)
		  {
		    double diff = (grad[i+i2][j+j2] - means[ii][i2+template_size][j2+template_size]);
		    res[i][j] += -diff * diff / vars[ii][i2+template_size][j2+template_size];
		  }
	    }

	unary_potentials[ii] = res;
      }
    return unary_potentials;
  }

public:
    SDoublePlane means[6], vars[6];
  int template_size;
};
