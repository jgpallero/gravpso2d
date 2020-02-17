%This script creates the mex files written in C
%Compiling using OpenMP: 0/1 -> no/yes
omp = 1;
%Check if I want OpenMP
if omp
    %Check if we are in GNU Octave or Matlab
    if exist('OCTAVE_VERSION')
        mex -O3 -fopenmp grav2d_GravityRectangle.c
        delete('grav2d_GravityRectangle.o');
    else
        %Check between MS Windows and no MS Windows
        if ispc~=0
            %Provides OpenMP support with MSVC options
            mex COMPFLAGS="$COMPFLAGS /O2 /openmp" grav2d_GravityRectangle.c
        else
            %Provides OpenMP support with GCC options
            mex CFLAGS="$CFLAGS -O3 -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"...
                grav2d_GravityRectangle.c
        end
    end
else
    %Check if we are in GNU Octave or Matlab
    if exist('OCTAVE_VERSION')
        %Serial version
        mex -O3 grav2d_GravityRectangle.c
        delete('grav2d_GravityRectangle.o');
    else
        %Check between MS Windows and no MS Windows
        if ispc~=0
            %Serial version
            mex COMPFLAGS="$COMPFLAGS /O2" grav2d_GravityRectangle.c
        else
            %Serial version
            mex CFLAGS="$CFLAGS -O3" grav2d_GravityRectangle.c
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (c) 2018, J.L.G. Pallero, jgpallero@gmail.com,
%                    J.L. Fernández Martínez, jlfm@uniovi.es
%                    Z. Fernández Muñiz, zulima@uniovi.es
%                    Sylvain Bonvalot, sylvain.bonvalot@ird.fr
%
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%- Redistributions of source code must retain the above copyright notice, this
%  list of conditions and the following disclaimer.
%- Redistributions in binary form must reproduce the above copyright notice,
%  this list of conditions and the following disclaimer in the documentation
%  and/or other materials provided with the distribution.
%- Neither the name of the copyright holders nor the names of its contributors
%  may be used to endorse or promote products derived from this software without
%  specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
%INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
%LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
%OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
%ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
