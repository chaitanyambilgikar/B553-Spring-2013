#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <map>


using namespace std;

// Definition for a 2D pixel point.
// Add additional fields if necessary.
// 
class Coordinate
{  
public:
  Coordinate() :row(0), col(0) {};
  double row, col;
};

// Ground truth nformation for training or testing image
class GTImage
{
public:
  string image_filename;
  Coordinate coords[6];
};

// Code for computing unary potential functions
//  Feel free to modify this if you want, but 
//  it should work as-is.
#include <UnaryPotential.h>


// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlay rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}





// Load ground truth coordinates and image names from file
//
vector<GTImage> load_groundtruth(const string &gt_filename)
{
  // read in training images, one at a time, and put the ground truth coordinates
  //  into a vector.
  ifstream ifs(gt_filename.c_str());
  
  vector<GTImage> training_images;
  while(ifs.good())
    {
      string image_number;
      ifs >> image_number;
      if(!ifs.good())	break;
      
      GTImage gt;
      gt.image_filename = string("imgs/image_") + image_number + ".png";
      for(int ii=0; ii<6; ii++)
	ifs >> gt.coords[ii].col >> gt.coords[ii].row;

      training_images.push_back(gt);
    }

  return training_images;
}



// Draw a given part configuration on an image, using simple
//  colored rectangles.
// Color code is: 
// left eye: bright green
// right eye: dark green
// nose: bright red
// left mouth: bright blue
// right mouth: dark blue
// chin: dark red
void draw_configuration(const char *output_fname, const char *sample_filename, const vector<Coordinate> &configuration)
{
  int g_code[] = {255, 128,   0,   0,   0, 0};
  int r_code[] = {0, 0,     255,   0,   0, 128};
  int b_code[] = {0, 0,       0, 255, 128, 0};

  SDoublePlane img = SImageIO::read_png_file(sample_filename);
  SDoublePlane r = img, g = img, b = img;
  for(int ii=0; ii<configuration.size(); ii++)
    {
      overlay_rectangle(g, configuration[ii].row-15, configuration[ii].col-15, configuration[ii].row+15, configuration[ii].col+15, g_code[ii], 3);
      overlay_rectangle(r, configuration[ii].row-15, configuration[ii].col-15, configuration[ii].row+15, configuration[ii].col+15, r_code[ii], 3);
      overlay_rectangle(b, configuration[ii].row-15, configuration[ii].col-15, configuration[ii].row+15, configuration[ii].col+15, b_code[ii], 3);
    }

  SImageIO::write_png_file(output_fname, r, g, b);
}




void learning(vector<GTImage> training_data, Coordinate *means[6], Coordinate *vars[6])
{
  
    //intialize all phi's to zero
    for (int k=0; k<6; k++)
    {
        for (int j=k+1; j<6; j++)         //means for combinations of parts, eg Mu: righteye,lefteye would be means[0][1], ie lower index always comes first
        {                            
          means[k][j].row=0;
          means[k][j].col=0;
          vars[k][j].row=0;
          vars[k][j].col=0;            //variances set up same way as means
          //mean_row[k][j]=0;
          //mean_col[k][j]=0;
          //var_row[k][j]=0;
          //var_col[k][j]=0;
        }
    }
    //sum up differences squared of possible combinations of parts to calculate means 
    for (int i=0; i<training_data.size(); i++)
    {
      for (int k=0; k<6; k++)
      {
        for (int j=k+1; j<6; j++)
        {
          means[k][j].row += (pow  ( (training_data[i].coords[k].row - training_data[i].coords[j].row) ,2)) ;
          means[k][j].col += (pow  ( (training_data[i].coords[k].col - training_data[i].coords[j].col) ,2)) ;       
        }
      }
    }
    //caclulate mean of differences by dividing sums by total number of elements
    for (int k=0; k<6; k++)
    {
    	for (int j=k+1; j<6; j++)
    	{
    		means[k][j].row = means[k][j].row/training_data.size();
    		means[k][j].col = means[k][j].col/training_data.size();
    	}
    }
    
    //calculations of sums to calculate variances
    for (int i=0; i<training_data.size(); i++)
    {
  	for (int k=0; k<6; k++)
  	{
  		for(int j=k+1; j<6; j++)
  		{
  			//vars[k][j].row += (pow( (pow  (training_data[i].coords[k].row - training_data[i].coords[j].row),2 ) ) - means[k][j].row) , 2) ;
  			//vars[k][j].col += (pow( (pow  ((training_data[i].coords[k].col - training_data[i].coords[j].col),2 ) ) - means[k][j].col) , 2) ;

        vars[k][j].row +=  pow(((pow((training_data[i].coords[k].row - training_data[i].coords[j].row), 2)) - means[k][j].row) , 2);
        vars[k][j].col +=  pow(((pow((training_data[i].coords[k].col- training_data[i].coords[j].col), 2)) - means[k][j].col) , 2);
  		}
  	}
    }
  
  //divide sums to get variances 
  for (int k=0; k<6; k++)
  {
  	for (int j=k+1; j<6; j++)
  	{
  		vars[k][j].row = vars[k][j].row/training_data.size();
  		vars[k][j].col = vars[k][j].col/training_data.size();
  	}
  }

  /*for (int k=0; k<6; k++)
  {
    for (int j=0; j<6; j++)
    {
      means[j][k].row = -means[k][j].row;
      means[j][k].col = -means[k][j].col;
      vars[j][k].row = -vars[k][j].row;
      vars[j][k].col = -vars[k][j].col;
      /*cout << "Mean row  of part[" << k << "][" << j << "] : " << means[k][j].row << endl;
      cout << "Mean col  of part[" << k << "][" << j << "] : " << means[k][j].col << endl;
      cout << "Variance row  of part[" << k << "][" << j << "] : " << vars[k][j].row << endl;
      cout << "Variance col  of part[" << k << "][" << j << "] : " << vars[k][j].col << endl;
    }
  }*/
}



SDoublePlane Phi_calculator(vector <vector<Coordinate> > &coor_map, Coordinate mean, Coordinate var, SDoublePlane my_potentials,int rows, int cols) //returns an array of coordinates, which is the best possible assignments for two parts i and j
{
  SDoublePlane phi = SDoublePlane(rows,cols);
  cout<<"\nReached here \n";
  for(int ix = 0; ix<rows; ix++)
  {
    for(int iy = 0; iy<cols; iy++)
    {
      double max_phi = -1e100;
      Coordinate max_coord;
      max_coord.row = 0;
      max_coord.col = 0;
      for (int jx = 0; jx<rows; jx++)
      {
        for (int jy = 0; jy<cols; jy++)
        {
          double t1 = pow(ix-jx-mean.row,2)/(2*var.row);
          double t2 = pow(iy-jy-mean.col,2)/(2*var.col);
          double temp = -t1-t2+my_potentials[jx][jy];
          //cout << "\n Phi  " << temp << endl;
          if(temp>max_phi)
          {
            max_phi = temp;
            max_coord.row = jx;
            max_coord.col = jy;
          }
        }
      }
      phi[ix][iy] = max_phi;
      coor_map[ix][iy].row = max_coord.row;
      coor_map[ix][iy].col = max_coord.col;
      //cout << "For parents co ords " << ix << " " << iy << " childs coords are " << max_coord.row << " " << max_coord.col << endl;
    }
  }
  return phi;

} 



void MAP_Inference (vector<SDoublePlane> my_potentials, string filename, Coordinate *means[6], Coordinate *vars[6]) 
{
//takes all training data loaded by load_ground_truth and means and vars meade by learning(). Also takes mypotentials which are already ready in main from where this method will be called.


	//we need unaries of leaves of graph B, which are right eye, nose, left mouth, right mouth, and chin. Their indexes are 1,2,3,4,5
	//we find max for each
	/*
  vector<Coordinate> best_configuration(6);
	double max_val[6];
	//initialize max_vals for each part with unary potentials and find and store max of those
	for(int ii=0; ii<6; ii++)
	{
		max_val[ii] = -1e100;
		for(int i=0; i<my_potentials[ii].rows(); i++)
			for(int j=0; j<my_potentials[ii].cols(); j++)
				if(my_potentials[ii][i][j] > max_val[ii])
				{
					best_configuration[ii].row = i;
					best_configuration[ii].col = j;
					max_val[ii] = my_potentials[ii][i][j];
					//cout << my_potentials[ii][i][j] <<endl;
					//cout << "max" << max_val << endl;
        }
  }
  */
	//getting details about the image used
  SDoublePlane image1 = SImageIO::read_png_file(filename.c_str());
  int rows = image1.rows();//number of rows in the image
  int cols = image1.cols();	//number of cols in the image

  
  //********************************************GRAPH B *********************************************************
  vector < vector <Coordinate> > coor_map_le_re(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_le_n(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_le_rm(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_le_lm(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_le_c(rows,vector<Coordinate>(cols));
  for (int ii = 0;ii<rows;ii++)
  {
    for(int jj = 0;jj<cols;jj++)
    {
      coor_map_le_re[ii][jj].row =0;
      coor_map_le_re[ii][jj].col = 0;

      coor_map_le_n[ii][jj].row =0;
      coor_map_le_n[ii][jj].col = 0;

      coor_map_le_rm[ii][jj].row =0;
      coor_map_le_rm[ii][jj].col = 0;

      coor_map_le_lm[ii][jj].row =0;
      coor_map_le_lm[ii][jj].col = 0;

      coor_map_le_c[ii][jj].row =0;
      coor_map_le_c[ii][jj].col = 0;
    }
  }
  cout<<"Mean of LE-RE "<<means[0][1].row<<"Vars of LE-RE "<<vars[0][1].row;



  SDoublePlane lefteye_righteye = Phi_calculator(coor_map_le_re, means[0][1], vars[0][1], my_potentials[1], rows, cols);
  SDoublePlane lefteye_nose = Phi_calculator(coor_map_le_n, means[0][2], vars[0][2], my_potentials[2], rows, cols);
  SDoublePlane lefteye_leftmouth = Phi_calculator(coor_map_le_lm, means[0][3], vars[0][3], my_potentials[3], rows, cols);
  SDoublePlane lefteye_rightmouth = Phi_calculator(coor_map_le_rm, means[0][4], vars[0][4], my_potentials[4], rows, cols);
  SDoublePlane lefteye_chin = Phi_calculator(coor_map_le_c, means[0][5], vars[0][5], my_potentials[5], rows, cols);

  SDoublePlane lefteye = SDoublePlane(rows, cols);
  vector<Coordinate> assignments_B(6);
  double temp = 0.0;
  double max = -1e100;
  for (int i=0; i<rows; i++)
  {
    for (int j=0; j<cols; j++)
    {
      temp = lefteye_righteye[i][j] + lefteye_nose[i][j] + lefteye_leftmouth[i][j] + lefteye_rightmouth[i][j] + lefteye_chin[i][j] + my_potentials[0][i][j];
      if (temp > max)
      {
        max = temp;
        assignments_B[0].row = i;
        assignments_B[0].col = j;
      }
    }
  } 

  cout << "left eye row and col " << assignments_B[0].row << " " << assignments_B[0].col << endl;

  assignments_B[1].row = coor_map_le_re[assignments_B[0].row][assignments_B[0].col].row;
  assignments_B[1].col = coor_map_le_re[assignments_B[0].row][assignments_B[0].col].col;

  cout << "right eye row and col " << assignments_B[1].row << " " << assignments_B[1].col << endl;

  assignments_B[2].row = coor_map_le_n[assignments_B[0].row][assignments_B[0].col].row;
  assignments_B[2].col = coor_map_le_n[assignments_B[0].row][assignments_B[0].col].col;

  cout << "nose row and col " << assignments_B[2].row << " " << assignments_B[2].col << endl;

  assignments_B[3].row = coor_map_le_lm[assignments_B[0].row][assignments_B[0].col].row;
  assignments_B[3].col = coor_map_le_lm[assignments_B[0].row][assignments_B[0].col].col;

  cout << "left mouth row and col " << assignments_B[3].row << " " << assignments_B[3].col << endl;

  assignments_B[4].row = coor_map_le_rm[assignments_B[0].row][assignments_B[0].col].row;
  assignments_B[4].col = coor_map_le_rm[assignments_B[0].row][assignments_B[0].col].col;

  cout << "right mouth row and col " <<  assignments_B[4].row << " " << assignments_B[4].col << endl;

  assignments_B[5].row = coor_map_le_c[assignments_B[0].row][assignments_B[0].col].row;
  assignments_B[5].col = coor_map_le_c[assignments_B[0].row][assignments_B[0].col].col;

  cout << "chin row and col " << assignments_B[5].row << " " << assignments_B[5].col << endl;

  draw_configuration("resultMARKOV1.png", filename.c_str(), assignments_B);  





  /*
  //***********************************************GRAPH C****************************************************
  vector < vector <Coordinate> > coor_map_n_le(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_n_re(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_n_lm(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_n_rm(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_n_c(rows,vector<Coordinate>(cols));

  for (int ii = 0;ii<rows;ii++)
  {
    for(int jj = 0;jj<cols;jj++)
    {
      coor_map_n_le[ii][jj].row =0;
      coor_map_n_le[ii][jj].col = 0;

      coor_map_n_re[ii][jj].row =0;
      coor_map_n_re[ii][jj].col = 0;

      coor_map_n_lm[ii][jj].row =0;
      coor_map_n_lm[ii][jj].col = 0;

      coor_map_n_rm[ii][jj].row =0;
      coor_map_n_rm[ii][jj].col = 0;

      coor_map_n_c[ii][jj].row =0;
      coor_map_n_c[ii][jj].col = 0;
    }
  }

  SDoublePlane nose_lefteye = Phi_calculator(coor_map_n_le, means[0][2], vars[0][2], my_potentials[0], rows, cols);
  SDoublePlane nose_righteye = Phi_calculator(coor_map_n_re, means[1][2], vars[1][2], my_potentials[1], rows, cols);
  SDoublePlane nose_leftmouth = Phi_calculator(coor_map_n_lm, means[2][3], vars[2][3], my_potentials[3], rows, cols);
  SDoublePlane nose_rightmouth = Phi_calculator(coor_map_n_rm, means[2][4], vars[2][4], my_potentials[4], rows, cols);
  SDoublePlane nose_chin = Phi_calculator(coor_map_n_c, means[2][5], vars[2][5], my_potentials[5], rows, cols);

  SDoublePlane nose = SDoublePlane(rows, cols);
  vector<Coordinate> assignments_B(6);
  double temp = 0.0;
  double max = -1e100;
  for (int i=0; i<rows; i++)
  {
    for (int j=0; j<cols; j++)
    {
      temp = nose_lefteye[i][j] + nose_righteye[i][j] + nose_leftmouth[i][j] + nose_rightmouth[i][j] + nose_chin[i][j] + my_potentials[2][i][j];
      if (temp > max)
      {
        max = temp;
        assignments_B[2].row = i;
        assignments_B[2].col = j;
      }
    }
  } 

  //cout << "left eye row and col " << assignments_B[0].row << " " << assignments_B[0].col << endl;

  assignments_B[0].row = coor_map_n_le[assignments_B[2].row][assignments_B[2].col].row;
  assignments_B[0].col = coor_map_n_le[assignments_B[2].row][assignments_B[2].col].col;

  //cout << "right eye row and col " << assignments_B[1].row << " " << assignments_B[1].col << endl;

  assignments_B[2].row = coor_map_n_re[assignments_B[2].row][assignments_B[2].col].row;
  assignments_B[2].col = coor_map_n_re[assignments_B[2].row][assignments_B[2].col].col;

  //cout << "nose row and col " << assignments_B[2].row << " " << assignments_B[2].col << endl;

  assignments_B[3].row = coor_map_n_lm[assignments_B[2].row][assignments_B[2].col].row;
  assignments_B[3].col = coor_map_n_lm[assignments_B[2].row][assignments_B[2].col].col;

  //cout << "left mouth row and col " << assignments_B[3].row << " " << assignments_B[3].col << endl;

  assignments_B[4].row = coor_map_n_rm[assignments_B[2].row][assignments_B[2].col].row;
  assignments_B[4].col = coor_map_n_rm[assignments_B[2].row][assignments_B[2].col].col;

  //cout << "right mouth row and col " <<  assignments_B[4].row << " " << assignments_B[4].col << endl;

  assignments_B[5].row = coor_map_n_c[assignments_B[2].row][assignments_B[2].col].row;
  assignments_B[5].col = coor_map_n_c[assignments_B[2].row][assignments_B[2].col].col;

  //cout << "chin row and col " << assignments_B[5].row << " " << assignments_B[5].col << endl;



  draw_configuration("resultMARKOV.png", filename.c_str(), assignments_B);  
   


  //************************************* GRAPH D ***************************************************
  vector < vector <Coordinate> > coor_map_re_le(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_rm_le(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_n_lm(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_n_rm(rows,vector<Coordinate>(cols));
  vector < vector <Coordinate> > coor_map_n_c(rows,vector<Coordinate>(cols));

  for (int ii = 0;ii<rows;ii++)
  {
    for(int jj = 0;jj<cols;jj++)
    {
      coor_map_n_le[ii][jj].row =0;
      coor_map_n_le[ii][jj].col = 0;

      coor_map_n_re[ii][jj].row =0;
      coor_map_n_re[ii][jj].col = 0;

      coor_map_n_lm[ii][jj].row =0;
      coor_map_n_lm[ii][jj].col = 0;

      coor_map_n_rm[ii][jj].row =0;
      coor_map_n_rm[ii][jj].col = 0;

      coor_map_n_c[ii][jj].row =0;
      coor_map_n_c[ii][jj].col = 0;
    }
  }
  




} 






/*
void accuracy(vector<GTImage> training_data, vector<GTImage> testing_answers and few more parameters)
{
  //MAP accuracy
  for i in testingdata.size
  {
    for (j=0; j<6; j++)
    {
      if ( abs(testingdata[i].coords[j].row - trainingdata[i].coords[j].row) < 20)
      {
        correctcounters ++;
      }
      if ( abs(testingdata[i].coords[j].col - trainingdata[i].coords[j].col) < 20)
      {
        correctcounters ++ //increase appropriate counters
      }
      total count = testingdata.size();
    }

  }

} */

}

//*************************************************************************
//This function calculates the phi(Li,Lj) for the loopy function
SDoublePlane Phi_calculator_loopy(/*vector <vector<Coordinate> > &coor_map,*/ Coordinate mean, Coordinate var,int rows, int cols) //returns an array of coordinates, which is the best possible assignments for two parts i and j
{
  SDoublePlane phi = SDoublePlane(rows,cols);
  cout<<"\nReached here \n";
  for(int ix = 0; ix<rows; ix++)
  {
    for(int iy = 0; iy<cols; iy++)
    {
      double max_phi = -1e100;
      //Coordinate max_coord;
      //max_coord.row = 0;
      //max_coord.col = 0;
      for (int jx = 0; jx<rows; jx++)
      {
        for (int jy = 0; jy<cols; jy++)
        {
          double t1 = pow(ix-jx-mean.row,2)/(2*var.row);
          double t2 = pow(iy-jy-mean.col,2)/(2*var.col);
          double temp = -t1-t2;
          if(temp>max_phi)
          {
            max_phi = temp;
            //max_coord.row = jx;
            //max_coord.col = jy;
          }
        }
      }
      phi[ix][iy] = max_phi;
      //coor_map[ix][iy].row = max_coord.row;
      //coor_map[ix][iy].col = max_coord.col;
    }
  }
  return phi;

} 
SDoublePlane add (SDoublePlane a,SDoublePlane b,int rows,int cols)
{
  SDoublePlane result = SDoublePlane(rows,cols);
  for (int i = 0; i < rows ; i++ )
  {
    for ( int j = 0; j < cols ; j++ )
    {
      result[i][j] += a[i][j] + b[i][j];
    }
  }
  return result;
}
void assign (SDoublePlane &lhs,SDoublePlane rhs,int rows,int cols)
{
  for (int i =0; i<rows;i++)
  {
    for (int j=0;j<cols;j++)
    {
      lhs[i][j] = rhs[i][j];
    }
  }
}

void loopy(vector <SDoublePlane>my_potentials, string filename, Coordinate *means[6], Coordinate *vars[6])
{
  /*
  Here the messages are from LE-LM, LM-CH, CH-N, CH-RM, N-RE, RM-RE,RE-LE
  */
  //getting details about the image used
  SDoublePlane image1 = SImageIO::read_png_file(filename.c_str());
  int rows = image1.rows();//number of rows in the image
  int cols = image1.cols(); //number of cols in the image
  //Initialize the messages
  /*
  The messages correspond to the model as follows:
  0 - Left Eye    -->   Left Mouth
  1 - Left Mouth  -->   Chin
  2 - Chin        -->   Nose
  3 - Chin        -->   Right Mouth
  4 - Nose        -->   Right Eye
  5 - Right Mouth -->   Right Eye
  6 - Right Eye   -->   Left Eye
  
  vector <SDoublePlane>messages (7,SDoublePlane(rows,cols));
  */
  vector <SDoublePlane>node_message (6,SDoublePlane(rows,cols));
  //give a default value to all the messages
  for ( int ii=0 ; ii<6 ; ii++ )
  {
    for ( int i=0 ; i<rows ; i++ )
    {
      for ( int j=0 ; j<cols ; j++ )
      {
        //if( ii != 6 )
        //{
          node_message[ii][i][j] = 1;
        //}
        //messages[ii][i][j] = 1;
      }
    }
  }
  //Main Loop
  for (int i = 0; i<4 ; i++)
  {
    for (int node = 0; node < 6; node++ )
   {
    cout<<"\n Computing messages sent by node "<<node;
    if ( node == 0 )
    {
      assign(node_message[3],add(node_message[3],add(add(node_message[1],my_potentials[0],rows,cols),Phi_calculator_loopy(means[0][3],vars[0][3],rows,cols),rows,cols),rows,cols),rows,cols);
      assign(node_message[1],add(node_message[1],add(add(node_message[3],my_potentials[0],rows,cols),Phi_calculator_loopy(means[0][1],vars[0][1],rows,cols),rows,cols),rows,cols),rows,cols);
      /*
      messages[0] = add(add(messages[6],my_potentials[0],rows,cols),Phi_calculator_loopy(means[0][3],vars[0][3],rows,cols),rows,cols);
      messages[6] = add(add(messages[0],my_potentials[0],rows,cols),Phi_calculator_loopy(means[0][1],vars[0][1],rows,cols),rows,cols);
      node_message[node] = add(messages[0],messages[6],rows,cols);
      */
    }
    if ( node == 1 )
    {
      assign(node_message[0] , add(add(node_message[0],node_message[2],rows,cols),add(add(node_message[4],my_potentials[1],rows,cols),Phi_calculator_loopy(means[0][1],vars[0][1],rows,cols),rows,cols),rows,cols),rows,cols);
      assign(node_message[2] , add(add(node_message[2],node_message[0],rows,cols),add(add(node_message[4],my_potentials[1],rows,cols),Phi_calculator_loopy(means[1][2],vars[1][2],rows,cols),rows,cols),rows,cols),rows,cols);
      assign(node_message[4] , add(add(node_message[4],node_message[0],rows,cols),add(add(node_message[2],my_potentials[1],rows,cols),Phi_calculator_loopy(means[1][4],vars[1][4],rows,cols),rows,cols),rows,cols),rows,cols);
    }
    if ( node == 2 )
    {
      assign(node_message[1] , add(add(node_message[1],node_message[5],rows,cols),add(my_potentials[2],Phi_calculator_loopy(means[1][2],vars[1][2],rows,cols),rows,cols),rows,cols),rows,cols);
      assign(node_message[5] , add(add(node_message[5],node_message[1],rows,cols),add(my_potentials[2],Phi_calculator_loopy(means[2][5],vars[2][5],rows,cols),rows,cols),rows,cols),rows,cols);
    }
    if ( node == 3 )
    {
      assign(node_message[0] , add(add(node_message[0],node_message[5],rows,cols),add(my_potentials[3],Phi_calculator_loopy(means[0][3],vars[0][3],rows,cols),rows,cols),rows,cols),rows,cols);
      assign(node_message[5] , add(add(node_message[5],node_message[0],rows,cols),add(my_potentials[3],Phi_calculator_loopy(means[3][5],vars[3][5],rows,cols),rows,cols),rows,cols),rows,cols);
    }
    if ( node == 4)
    {
      assign(node_message[1] , add(add(node_message[1],node_message[5],rows,cols),add(my_potentials[4],Phi_calculator_loopy(means[1][4],vars[1][4],rows,cols),rows,cols),rows,cols),rows,cols);
      assign(node_message[5] , add(add(node_message[5],node_message[1],rows,cols),add(my_potentials[4],Phi_calculator_loopy(means[4][5],vars[4][5],rows,cols),rows,cols),rows,cols),rows,cols);
    }
    if( node == 5 )
    {
      assign(node_message[2] , add(add(node_message[2],node_message[3],rows,cols),add(add(node_message[4],my_potentials[5],rows,cols),Phi_calculator_loopy(means[2][5],vars[2][5],rows,cols),rows,cols),rows,cols),rows,cols);
      assign(node_message[3] , add(add(node_message[3],node_message[2],rows,cols),add(add(node_message[4],my_potentials[5],rows,cols),Phi_calculator_loopy(means[3][5],vars[3][5],rows,cols),rows,cols),rows,cols),rows,cols);
      assign(node_message[4] , add(add(node_message[4],node_message[2],rows,cols),add(add(node_message[3],my_potentials[5],rows,cols),Phi_calculator_loopy(means[4][5],vars[4][5],rows,cols),rows,cols),rows,cols),rows,cols);
    }
  }
  }
  
  vector <Coordinate>best_configuration(6);
  for ( int node = 0; node < 6 ; node++ )
  {
    double temp = -1e10000;
    for (int i = 0; i<rows ; i++)
    {
      for ( int j = 0; j<cols ; j++ )
      {
        //cout<<node_message[node][i][j]<<" ";
        if (node_message[node][i][j] > temp)
        {
          temp = node_message[node][i][j];
          best_configuration[node].row = i;
          best_configuration[node].col = j;
        }
      }
      //cout<<"\n";
      //cout<<"\nBest Config for node "<<node<<" is "<<best_configuration[node].row<<","<<best_configuration[node].col;
    }
    //cout<<"\n\n\n";
  }
  draw_configuration("resultLOOPY.png", filename.c_str(), best_configuration);
  //cout<<"\a\a\a\a\a\a";

}


int main(int argc, char *argv[])
{
  try {
    if(!(argc == 3)) {
      cerr << "usage: " << argv[0] << " training_gt_file testing_gt_file" << endl;
      return 1;
    }

    string train_gt_filename(argv[1]), test_gt_filename(argv[2]);

    // This is a parameter used to compute the unary potentials. It specifies
    //  how big the part appearance models should be.
    const int template_size=15;

    // Load in the training data
    vector<GTImage> training_data = load_groundtruth(train_gt_filename);


    Coordinate *means[6]; //parameters for training
    Coordinate *vars[6];
    for (int i=0; i<6; i++)
    {
      means[i] = new Coordinate[6];
      vars[i] = new Coordinate[6];
    }
    learning(training_data, means, vars);



    // Now set up the unary potential functions
    UnaryPotential unary_potentials(training_data, template_size);
    cerr << "Done learning unary potentials" << endl << endl;

    
    

    // Now run naive inference on one of the sample images
    //  
    string sample_filename = "imgs/image_0001.png";
    cerr << "Running naive inference on " << sample_filename << "..." << endl;

    // Compute unary potentials for this image
    vector<SDoublePlane> my_potentials = unary_potentials.compute_unary_logpotential(sample_filename);

    //MAP_Inference (my_potentials, "imgs/image_0001.png", means, vars);

    loopy(my_potentials, "imgs/image_0001.png", means, vars);


    // Naive inference just involves choosing best coordinate in each individual
    //  unary potential function.
    const int part_count=6;
    vector<Coordinate> best_configuration(part_count);
    for(int ii=0; ii<part_count; ii++)
      {
	double max_val=-1e100;
	for(int i=0; i<my_potentials[ii].rows(); i++)
	  for(int j=0; j<my_potentials[ii].cols(); j++)
	    if(my_potentials[ii][i][j] > max_val)
	      {
		best_configuration[ii].row = i;
		best_configuration[ii].col = j;
		max_val = my_potentials[ii][i][j];
    //cout << my_potentials[ii][i][j] <<endl;
    //cout << "max" << max_val << endl;
	      }
      }

    // Visualize result by drawing detected parts on the image, and 
    //  outputing to a new image file.
    //
    draw_configuration("naive_result.png", sample_filename.c_str(), best_configuration);

  } catch(std::string &str) {
    cerr << str << endl;
  }

  return 0;
}