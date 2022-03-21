function pos = getMotion( t, varargin )
% getMotion Simulate 3D motion for use in MRI PSF calculation
%
% pos = getMotion( t, 'param1', value1, 'param2', value2, ...)
%
%           returns Cartesian coordinates of a moving point
%           at all time bins in t.
%
% Input
% =====
%   t       - (mandatory) vector of time points, unit is [s]
%
% Parameters to control the modelled motion can be given as key/value pairs.
%
%
% The actual motion model has to be altered inside this function!
%
%
% Output
% ======
%   pos     -   position (x,y,z) for each point in time
%               size(pos) = [ length(t)   3 ]
%
% Example
% =======
%
%   Simulate a sinusoidal motion in the x-direction (== readout direction)
%   with 8.5 Px amplitude and a frequency of 5 Hz using 1000 time bins
%   in the range 0...5 seconds:
%
%       t = linspace( 0, 5, 1000 );
%       pos = getMotion( t, 'sinXampl', 8.5, 'sinXfreq', 5 );
%
%   Valid parameters are any out of:
%
%   'constX'  , 'constY'  , 'constZ'    for a constant offset
%   'shiverX' , 'shiverY' , 'shiverZ'   for random motion (shiverX/Y/Z defining the std deviation)
%   'sinXampl', 'sinYampl', 'sinZampl'  Amplitudes of sinusoidal motion on either axis
%   'sinXfreq', 'sinYfreq', 'sinZfreq'    Frequencies accordingly
%   'cosXampl', 'cosYampl', 'cosZampl'  Amplitudes of cosinusoidal motion on either axis
%   'cosXfreq', 'cosYfreq', 'cosZfreq'    Frequencies accordingly
%   'respiXampl', 'respiYampl', 'respiZampl'  Amplitudes of respiratory motion on either axis
%   'respirFreq'                                Respiration frequency (default: 10/60s)
%   'heartXampl', 'heartYampl', 'heartZampl'  Amplitudes of cardiac motion on either axis
%   'heartFreq'                                 Cardiac frequency (default: 75/60s)
%
%   To plot the outcome of the simulation:
%
%   plot( t, pos, 'LineWidth', 2), legend('\Delta X', '\Delta Y', '\Delta Z')

% michael.voelker@mr-bavaria.de, 2013

    %  parse input as key-value pairs
    %
    p = inputParser;
    isRealScalar = @(x) validateattributes( x, {'numeric'}, {'nonempty', 'scalar', 'real', 'finite'} );
    % addParamValue( 'varName'    , default, checkFunction   );
    p.addParamValue( 'constX'     ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'constY'     ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'constZ'     ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'shiverX'    ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'shiverY'    ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'shiverZ'    ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'sinXampl'   ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'sinYampl'   ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'sinZampl'   ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'sinXfreq'   ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'sinYfreq'   ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'sinZfreq'   ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'cosXampl'   ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'cosYampl'   ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'cosZampl'   ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'cosXfreq'   ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'cosYfreq'   ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'cosZfreq'   ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'respiXampl' ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'respiYampl' ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'respiZampl' ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'respirFreq' ,       0, @(x) isRealScalar(x)  ); %  mnh set to 0
%     p.addParamValue( 'respirFreq' ,     0.1, @(x) isRealScalar(x)  ); %  default: 0.1 Hz = 6 / min
    p.addParamValue( 'heartXampl' ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'heartYampl' ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'heartZampl' ,       0, @(x) isRealScalar(x)  );
    p.addParamValue( 'heartFreq'  ,       0, @(x) isRealScalar(x)  );%  mnh set to 0
%     p.addParamValue( 'heartFreq'  ,   75/60, @(x) isRealScalar(x)  ); %  default: 75 / min

    % Save a list of all available parameters which is easy to reference
    availParam = '';
    for n = 1:length(p.Parameters)
        availParam = [availParam  p.Parameters{n} '  '];
        if mod(n,6) == 0
            availParam = [ availParam  sprintf('\n')];
        end
    end

    % Generate useful info how to use getMotion() and which parameters are currently implemented
    helpStr = sprintf(['Try e.g.\n'    ...
                    '\t getMotion( t, ''sinXampl'', 8.5, ''sinXfreq'', 5 )\n'   ...
                    'to simulate sinusoidal motion in x-direction\n'    ...
                    'with 8.5 Px amplitude and 5 Hz frequency.\n\n'   ...
                    'The currently available parameters are:\n\n'     ...
                    '%s\n\n'     ...
                    'For more details:\n'   ...
                    '\thelp getMotion\n' ], availParam );
    switch nargin
        case 0
            error( 'getMotion:NoInput', [ 'First input: vector of time points in seconds.\n'    ...
                                           'further inputs: key/value pairs to specify the motion to be modelled.\n\n'  ...
                                           '%s' ], helpStr )
        case 1
            warning('getMotion:NoModel', ['No parameters were given to specify the motion to be modelled.\n'    ...
                    '%s' ], helpStr )
        otherwise
            % Well, things seem to be fine, then.
    end

    validateattributes( t, {'single', 'double'}, {'nonempty', 'vector', 'real', 'finite', 'nonnegative'} )

    t = t(:);       % make t a column vector
    Nt = numel( t );

    try
        p.parse( varargin{:} )
    catch err
        [a, b] = regexp( err.message, 'did not match any valid parameter' );
        if isempty(a) || isempty(b)
            additionalInfo = '';
        else
            additionalInfo = sprintf( 'You might have a simple typo in the parameter name.\n\n=== Use any name from this list (not case sensitive) ===\n%s\n', availParam );
        end
        error( 'getMotion:Input', '\n%s\n%s', err.message, additionalInfo )
    end
    r = p.Results;         % save the results

    model1_X = 0;   model2_X = 0;   model3_X = 0;   % \
    model1_Y = 0;   model2_Y = 0;   model3_Y = 0;   %   just be sure they exist
    model1_Z = 0;   model2_Z = 0;   model3_Z = 0;   % /


    %
    %   ACTUAL MOTION MODELS ARE HERE
    %

    
    % ( === random displacements ==========================================
    %
    shiverX  = r.shiverX * randn( Nt, 1);
    shiverY  = r.shiverY * randn( Nt, 1);
    shiverZ  = r.shiverZ * randn( Nt, 1);
    % ) ===================================================================


    % ( === simple periodical displacements ===============================
    %
    periodicalX =   r.sinXampl * sin( 2*pi * r.sinXfreq * t )       ...
                  + r.cosXampl * cos( 2*pi * r.cosXfreq * t );

    periodicalY =   r.sinYampl * sin( 2*pi * r.sinYfreq * t )       ...
                  + r.cosYampl * cos( 2*pi * r.cosYfreq * t );

    periodicalZ =   r.sinZampl * sin( 2*pi * r.sinZfreq * t )       ...
                  + r.cosZampl * cos( 2*pi * r.cosZfreq * t );
    % ) ===================================================================



    % ( === example for an explicit motion model  =========================
    % Model for respiration and heart beat
    % --> http://dx.doi.org/10.1002/mrm.22031
    %
        respiX = r.respiXampl * cos( 2*pi * r.respirFreq * t ).^(2*3);
        respiY = r.respiYampl * cos( 2*pi * r.respirFreq * t ).^(2*3);
        respiZ = r.respiZampl * cos( 2*pi * r.respirFreq * t ).^(2*3);

        heartX = r.heartXampl * cos( 2*pi * r.heartFreq * t ).^(2*2);
        heartY = r.heartYampl * cos( 2*pi * r.heartFreq * t ).^(2*2);
        heartZ = r.heartZampl * cos( 2*pi * r.heartFreq * t ).^(2*2);

        model1_X =  heartX  -  respiX;
        model1_Y =  heartY  -  respiY;
        model1_Z =  heartZ  -  respiZ;

        clear respi heart
    % ) ===================================================================



    % ( === Your model(s) here!  ==========================================
    %
    %
       % model2_X =      ;
       % model2_Y =      ;
       % model2_Z =      ;

       % model3_X =      ;
       % model3_Y =      ;
       % model3_Z =      ;
    % ) ===================================================================



    % ( === combine all kinds of motion ========================================================
    %
        deltaX  =  r.constX  +  shiverX  +  periodicalX  +  model1_X   +  model2_X   +  model3_X;
        deltaY  =  r.constY  +  shiverY  +  periodicalY  +  model1_Y   +  model2_Y   +  model3_Y;
        deltaZ  =  r.constZ  +  shiverZ  +  periodicalZ  +  model1_Z   +  model2_Z   +  model3_Z;
    % ) ========================================================================================


    % store in a common array with single precision to save memory:
    pos = single( [ deltaX, deltaY, deltaZ ] );

end % of getMotion()
