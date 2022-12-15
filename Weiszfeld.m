%Copyright (c) 2017, Francesco
%All rights reserved.

%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%* Redistributions of source code must retain the above copyright notice, this
%  list of conditions and the following disclaimer.
%
%* Redistributions in binary form must reproduce the above copyright notice,
%  this list of conditions and the following disclaimer in the documentation
%  and/or other materials provided with the distribution
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
%FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
%SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
%CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
%OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function output_structure = Weiszfeld(input_structure)
%
% This function numerically calculates the geometric mean of a
% N-Dimensional set of points using the Wieszfeld's algorithm
%
% INPUT: Structure
% -- [REQUIRED] input_strucuture.Data = Input data matrix. Each row is a point, each
% column a dimension
% -- input_structure.RelTol = Relative tolerance for stopping the search.
% Default: 0.001
% -- input_structure.x0 = A vector with the initial point. If not provided,
% it is automatically calcualted based on the centroid of the original
% series of points


%% Default values
RelTolDefault = 0.001 ;
expectedIterations = 20 ;

%% Reading inputs
% Reading the data matrix 
data = input_structure.Data ;
% Check if the RelTol field is provided. Otherwise, use the default value
if any(strcmp(fields(input_structure),'RelTol'))
    relTol = input_structure.RelTol ;
else
    relTol = RelTolDefault ;
end
% Check if a starting point is provided. Otherwise, calculate it
if any(strcmp(fields(input_structure),'x0'))
    x0 = input_structure.x0 ;
else
    x0 = mean(data,1) ;
end

%% Calculating some useful parameters
[nPoints, nDimensions] = size(data) ;

% Initialize the relative difference
eps = 1 ;
counter = 0 ;
% Initialize the matrix storing all iterations. 
xTemp = NaN([expectedIterations , nDimensions]) ;
xTemp(1,:) = x0 ;
%% Iterations
while eps > relTol
    counter = counter + 1 ;
    weights = sum((data - xTemp(counter,:)).^2, 2).^-0.5 ;
    temp = sum(data .* (weights .* ones(nPoints, nDimensions)), 1) / sum(weights) ;
    xTemp(counter+1, :) = temp ;
    eps = (sum((xTemp(counter+1,:) - xTemp(counter,:)).^2))^0.5 ;
end

%% Post compute
% Compute the difference at the last computation
err = sum(sum((data - xTemp(counter+1,:)).^2,2).^0.5) / nPoints;

output_structure.xMedian = xTemp(counter+1,:) ;
output_structure.err = err ;
output_strucutre.tol = eps ;

end
