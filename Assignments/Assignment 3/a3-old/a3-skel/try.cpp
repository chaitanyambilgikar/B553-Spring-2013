  for (int ix =0;ix<rows;ix++)
  {
    for(int iy=0;iy<cols;iy++)
    {
      double max_phi = -1e100;
      Coordinate max_cord;
      for (int jx=0;jx<rows;jx++)
      {
        for (int jy=0;jy<cols;jy++)
        {
          double temp = (-( (pow(ix-jx-mean.row,2))/(2*var.row) ) - ( (pow(iy-jy-mean.col,2))/(2*var.col) )) + my_potentials[jx][jy];
          if(temp>max_phi)
          {
            max_phi = temp;
            max_cord.row = jx;
            max_cord.col = jy;
          }
        }
      }
      phi[ix][iy] = max_phi;
      coor_map[ix][iy].row = max_cord.row;
      coor_map[ix][iy].col = max_cord.col;

    }
  }