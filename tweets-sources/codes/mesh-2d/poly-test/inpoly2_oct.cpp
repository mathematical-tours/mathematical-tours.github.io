
// a fast pre-sorted variant of the crossing-number test for 
// INPOLY2.m 

//----------------------------------------------------------
//  Darren Engwirda : 2017 --
//  Email           : de2363@columbia.edu
//  Last updated    : 17/07/2017
//----------------------------------------------------------

#include <octave/oct.h>

DEFUN_DLD (inpoly2_oct, args, ,
  "-- INPOLY2-OCT: low-level routine called by INPOLY2 to \n"
  "-- compute point-in-polygon queries.")
{    
    octave_value_list rval;
    
    int const nargin = args.length () ;
    if (nargin != +3)
    {
        print_usage () ;
        return rval ;
    }
    
    Matrix const vert (
        args(0).matrix_value ()) ;
    Matrix const node (
        args(1).matrix_value ()) ;
    Matrix const edge (
        args(2).matrix_value ()) ;
    
    if (error_state) return rval ;

    octave_idx_type const nvrt 
        = vert.rows () ;
    octave_idx_type const nnod 
        = node.rows () ;
    octave_idx_type const nedg 
        = edge.rows () ;

//---------------------------------- init. crossing no. bool
    boolMatrix stat(nvrt, 1) ;
    
    octave_idx_type vpos ;
    for (vpos = +0; vpos != nvrt; ++vpos)
    {
        stat(vpos) = false ;
    }

//---------------------------------- loop over polygon edges
    octave_idx_type epos ;
    for (epos = +0; epos != nedg; ++epos)
    {
        octave_idx_type const inod 
            = edge(epos,0) - 1 ;
        octave_idx_type const jnod 
            = edge(epos,1) - 1 ;

    //------------------------------ calc. edge bounding-box
        double yone = node (inod,1) ;
        double ytwo = node (jnod,1) ;
        double xone = node (inod,0) ;
        double xtwo = node (jnod,0) ;
        
        double ydel = ytwo - yone ;
        double xdel = xtwo - xone ;

        double xmin = xone < xtwo 
                    ? xone : xtwo ;
                    
    //------------------------------ find VERT(IPOS,2)<=YONE
        octave_idx_type ilow = +0 ; 
        octave_idx_type iupp = nvrt-1;
        octave_idx_type imid = +0 ;
        
        while (ilow < iupp - 1)
        {
            imid = ilow 
                + (iupp - ilow) / 2;
            
            if (vert(imid,1) < yone)
            {
                ilow = imid ;
            }
            else
            {
                iupp = imid ;
            }           
        }
        
    //------------------------------ calc. edge-intersection
        octave_idx_type vpos = ilow+1 ;
        for ( ; vpos != nvrt; ++vpos)
        {
            double xpos = vert(vpos,0);
            double ypos = vert(vpos,1);
            
            if (ypos <  ytwo)
            {
                if (xpos >= xmin)
                {
                    if (
                    ydel* (xpos-xone) <
                    xdel* (ypos-yone) )
                    {
                    stat(vpos) =
                       ! stat(vpos) ;
                    }
                }
                else
                {
                    stat(vpos) =
                       ! stat(vpos) ;
                }
            }
            else
            {
                break ;           // done -- due to the sort
            }
        }

    }

    rval(0) = stat;
    
    return rval;
}



