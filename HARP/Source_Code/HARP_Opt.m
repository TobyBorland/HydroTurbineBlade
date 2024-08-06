%=============================== HARP_Opt ================================%
% This is the main executable file, which initializes the GUI, initializes
% the genetic algorithm, and performs other various functions. 
%=========================================================================%

%=========================================================================%
% Initializes the GUI, DO NOT EDIT THIS SECTION
function varargout = HARP_Opt(varargin)

% HARP_Opt M-file for HARP_Opt.fig
%      HARP_Opt, by itself, creates a new HARP_Opt or raises the existing
%      singleton*.
%
%      H = HARP_Opt returns the handle to a new HARP_Opt or the handle to
%      the existing singleton*.
%
%      HARP_Opt('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HARP_Opt.M with the given input arguments.
%
%      HARP_Opt('Property','Value',...) creates a new HARP_Opt or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HARP_Opt_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HARP_Opt_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help HARP_Opt
% Last Modified by GUIDE v2.5 19-Jun-2010 15:18:22
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HARP_Opt_OpeningFcn, ...
                   'gui_OutputFcn',  @HARP_Opt_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
%=========================================================================%

%=========================================================================%
% Executes just before HARP_Opt is made visible.
function HARP_Opt_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HARP_Opt (see VARARGIN)

% Choose default command line output for HARP_Opt
handles.output = hObject;

%load the Background Image into the GUI
LogoImage = imread('Input_Files\Templates\NREL_Logo_Small_Blue.jpg');
axes(handles.GUI_Image);
image(LogoImage); axis off;
%set background color
set(gcf,'color',[0.933 0.933 0.933]);

% Update handles structure for the UI Button Groups
set(handles.RotSpdCtrl,'SelectionChangeFcn',@RotSpdCtrl_SelectionChangeFcn);
set(handles.BldPitCtrl,'SelectionChangeFcn',@BldPitCtrl_SelectionChangeFcn);
set(handles.OptimMethod,'SelectionChangeFcn',@OptimMethod_SelectionChangeFcn);
set(handles.ProbType,'SelectionChangeFcn',@ProbType_SelectionChangeFcn);
set(handles.TipLoss,'SelectionChangeFcn',@TipLoss_SelectionChangeFcn);
set(handles.HubLoss,'SelectionChangeFcn',@HubLoss_SelectionChangeFcn);
set(handles.Swirl,'SelectionChangeFcn',@Swirl_SelectionChangeFcn);
set(handles.AdvBrake,'SelectionChangeFcn',@AdvBrake_SelectionChangeFcn);
set(handles.IndProp,'SelectionChangeFcn',@IndProp_SelectionChangeFcn);
set(handles.AIDrag,'SelectionChangeFcn',@AIDrag_SelectionChangeFcn);
set(handles.TIDrag,'SelectionChangeFcn',@TIDrag_SelectionChangeFcn);
set(handles.DecST,'SelectionChangeFcn',@DecST_SelectionChangeFcn);
set(handles.Correct_3D,'SelectionChangeFcn',@Correct_3D_SelectionChangeFcn);
set(handles.ThickDist,'SelectionChangeFcn',@ThickDist_SelectionChangeFcn);
set(handles.ElementSpacing,'SelectionChangeFcn',@ElementSpacing_SelectionChangeFcn);
guidata(hObject, handles);
initialize_gui(hObject, handles, false);
%=========================================================================%

%=========================================================================%
% Shouldn't ever need to edit the function "varargout"
% --- Outputs from this function are returned to the command line.
function varargout = HARP_Opt_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%=========================================================================%

%=========================================================================%
function initialize_gui(fig_handle, handles, isreset)

%this section allows you to set the default values for the GUI inputs

%Turbine Configuration
    handles.SpdCtrl='1';
set(handles.RotSpdCtrl, 'SelectedObject', handles.SpdCtrl1);
    handles.PitCtrl='1';
set(handles.BldPitCtrl, 'SelectedObject', handles.PitCtrl1);
set(handles.NumBlade,'String','3');
set(handles.NumSeg,'String','30');
set(handles.RatedPwr,'String','1000');
set(handles.RotorDia,'String','50');
set(handles.HubDia,'String','2');
set(handles.HubHt,'String','75');
set(handles.PreCone,'String','0');
set(handles.ShaftTilt,'String','0');
set(handles.OmgMin,'String','4');
set(handles.OmgMax,'String','32.5');
    RPM = str2num(get(handles.OmgMax,'String'));
    R = 0.5*str2num(get(handles.RotorDia,'String'));
set(handles.TipSpd,'String',num2str(RPM*R*pi/30,'%6.1f'));

%Fluid Properties
set(handles.Rho,'String','1.225');
set(handles.KinVisc,'String','1.464e-5');
set(handles.SpdSt,'String','2');
set(handles.SpdEnd,'String','26');
set(handles.SpdDel,'String','0.5');

%Hydrokinetic Turbine Cavitation Inputs
set(handles.CheckCavit,'Value',0)
    if get(handles.CheckCavit,'Value')
        set(handles.VapPres,'Enable','on');
        set(handles.WatDepth,'Enable','on');
        set(handles.Patm,'Enable','on');
        set(handles.CavSF,'Enable','on');
     else
        % Checkbox is not checked-take approriate action
        set(handles.VapPres,'Enable','off');
        set(handles.WatDepth,'Enable','off');
        set(handles.Patm,'Enable','off'); 
        set(handles.CavSF,'Enable','off'); 
    end
set(handles.VapPres,'String','2500');
set(handles.WatDepth,'String','10');
set(handles.Patm,'String','101325');
set(handles.CavSF,'String','1.0');

%Optimization Objective
set(handles.OptimMethod,'SelectedObject',handles.OptAEP);
    handles.OptMethod = '1';
set(handles.ProbType, 'SelectedObject', handles.ProbDist1);
    handles.ProbDist='1';
set(handles.ProbDist1,'Enable','on');
set(handles.ProbDist2,'Enable','on');
set(handles.ProbDist3,'Enable','on');
set(handles.FlowFile_Button,'Enable','off');
set(handles.ProbDist4,'Enable','off');
set(handles.U_box,'Enable','on');
set(handles.k_box,'Enable','off');
set(handles.c_box,'Enable','off');
set(handles.Weib_U_box,'Enable','off');
set(handles.U_box,'String','6.03');
set(handles.k_box,'String','1.91');
set(handles.c_box,'String','6.8');
    cc = str2num(get(handles.c_box,'String'));
    kk = str2num(get(handles.k_box,'String'));
    Weib_Umean = cc*gamma(1+1/kk);
set(handles.Weib_U_box,'String',num2str(Weib_Umean,'%6.2f'));

%Structural Optimization inputs
set(handles.Einput,'String','27.6');
set(handles.allowableStrain,'String','3000');
set(handles.matDensity,'String','1800');
set(handles.STmin,'String','1');
set(handles.STdel,'String','0.2');
set(handles.SFstruct,'String','1.2');
set(handles.StructuralOpt,'Value',0)
set(handles.Einput,'Enable','off');
set(handles.allowableStrain,'Enable','off');
set(handles.matDensity,'Enable','off');
set(handles.STmin,'Enable','off');
set(handles.STdel,'Enable','off');
set(handles.SFstruct,'Enable','off');
set(handles.True12,'Enable','off');
set(handles.DecST, 'SelectedObject', handles.True12);
    handles.DecreasingST = '1';
set(handles.False12,'Enable','off');
set(handles.ParetoFraction,'Enable','off');
set(handles.EliteCount,'Enable','on');
set(handles.RecFail,'Value',0);

%blade element spacing
set(handles.ElementSpacing,'SelectedObject', handles.equal_spacing);
    handles.Elm_spacing = '0';

%Model Circular Root
set(handles.CircularRoot,'Value',0)
    if get(handles.CircularRoot,'Value')
        set(handles.minRtChord,'Enable','on');
        set(handles.maxRtChord,'Enable','on');
        set(handles.RtTranSt,'Enable','on');
        set(handles.RtTranEnd,'Enable','on');
    else
        % Checkbox is not checked-take approriate action
        set(handles.minRtChord,'Enable','off');
        set(handles.maxRtChord,'Enable','off');
        set(handles.RtTranSt,'Enable','off');
        set(handles.RtTranEnd,'Enable','off');
    end
set(handles.minRtChord,'String','0.8');
set(handles.maxRtChord,'String','1.5');
set(handles.RtTranSt,'String',num2str(str2num(get(handles.HubDia,'String'))/str2num(get(handles.RotorDia,'String'))+0.01));
set(handles.RtTranEnd,'String','0.25');

%Define the radial locations of the control points: using half-cosine_spacing spacing
   if get(handles.CircularRoot,'Value') == 0
   rad_CP = hcosspace(0.5*str2num(get(handles.HubDia,'String')),0.5*str2num(get(handles.RotorDia,'String')),5,0);    
   else
   rad_CP = hcosspace(str2num(get(handles.RtTranEnd,'String'))*0.5*str2num(get(handles.RotorDia,'String')),0.5*str2num(get(handles.RotorDia,'String')),5,0);      
   end
set(handles.radCP1,'String',num2str(rad_CP(1),4));
set(handles.radCP2,'String',num2str(rad_CP(2),4));
set(handles.radCP3,'String',num2str(rad_CP(3),4));
set(handles.radCP4,'String',num2str(rad_CP(4),4));
set(handles.radCP5,'String',num2str(rad_CP(5),4));

%BEM Algorithm Configuration
    handles.tiploss = 'True';
set(handles.TipLoss, 'SelectedObject', handles.True4);
    handles.hubloss = 'True';
set(handles.HubLoss, 'SelectedObject', handles.True5);
    handles.swirl = 'True';
set(handles.Swirl, 'SelectedObject', handles.True6);
    handles.advbrake = 'True';
set(handles.AdvBrake, 'SelectedObject', handles.True8);
    handles.indprop = 'True';
set(handles.IndProp, 'SelectedObject', handles.True9);
    handles.aidrag = 'True';
set(handles.AIDrag, 'SelectedObject', handles.True10);
    handles.tidrag = 'True';
set(handles.TIDrag, 'SelectedObject', handles.True11);
set(handles.MaxIter,'String','1000');
set(handles.ATol,'String','1.0e-6');

%Stall Delay Models
set(handles.Correct_3D, 'SelectedObject', handles.Correct_3D_0);
    handles.dragmodel = '0';
set(handles.Correct_3D_1,'Enable','off');
set(handles.Correct_3D_2,'Enable','off');

%Thickness distribution
set(handles.ThickDist, 'SelectedObject', handles.ThickMethod2);
    handles.ThickMethod = '2';

%Optimization Configuration
set(handles.SeedInitPop,'Value',0);
set(handles.NumGens,'String','150');
set(handles.PopSize,'String','200');
set(handles.EliteCount,'String','1');
set(handles.ParetoFraction,'String','0.5');
set(handles.ParetoFraction,'Enable','off');
set(handles.CrossFrc,'String','0.25');
set(handles.GATol,'String','1.0e-6');
set(handles.TwLB1,'String','-10');set(handles.TwLB2,'String','-10');set(handles.TwLB3,'String','-10');set(handles.TwLB4,'String','-10');set(handles.TwLB5,'String','-10');
set(handles.TwUB1,'String','40');set(handles.TwUB2,'String','40');set(handles.TwUB3,'String','40');set(handles.TwUB4,'String','40');set(handles.TwUB5,'String','40');
set(handles.CLB1,'String','0.5');set(handles.CLB2,'String','0.1');set(handles.CLB3,'String','0.1');set(handles.CLB4,'String','0.1');set(handles.CLB5,'String','0.1');
set(handles.CUB1,'String','1.5');set(handles.CUB2,'String','1.5');set(handles.CUB3,'String','1');set(handles.CUB4,'String','1');set(handles.CUB5,'String','0.25');
set(handles.ThickLB,'String','0.25 0.5 NaN');
set(handles.ThickUB,'String','0.35 NaN NaN');
set(handles.Family,'String','FFA');
set(handles.StepVals,'String','30 24 21');
set(handles.Filename,'String','Optimization_Run1');

UpdateEstimates(handles) %calls the function to estimate AEP and CF
set(handles.minD,'String','25');
set(handles.maxD,'String','75');
set(handles.minP,'String','500');
set(handles.maxP,'String','1500');

% Update handles structure
guidata(handles.figure1, handles);
%=========================================================================%

%=========================================================================%
%This code writes the template WT_Perf files: Input Configuration.wtp
%                                             Aerodynamic Inputs.wtp
%                                             Output Configuration.wtp

% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Open the files for writing, flushing any previous data in the files
fid1=fopen('Input_Files\Templates\Input_Configuration.wtp','wt');
fid2=fopen('Input_Files\Templates\Aerodynamic_Inputs.wtp','wt');
fid3=fopen('Input_Files\Templates\Output_Configuration.wtp','wt');

%Build the 'Input Configuration.wtp' text file
fprintf(fid1,'-----  WT_Perf Input File  -----------------------------------------------------\n');
fprintf(fid1,'This file (Non-optimal Solution) generated automatically by HARP_Opt during intermediate iterations \n');
fprintf(fid1,'Compatible with WT_Perf v3.00f                                                  \n');
fprintf(fid1,'-----  Input Configuration  ----------------------------------------------------\n');
fprintf(fid1,'%-21s','False');
fprintf(fid1,'Echo:                      Echo input parameters to "<rootname>.ech"?\n');
fprintf(fid1,'%-21s','True ');
fprintf(fid1,'DimenInp:                  Turbine parameters are dimensional?\n');
fprintf(fid1,'%-21s','True ');
fprintf(fid1,'Metric:                    Turbine parameters are Metric (MKS vs FPS)?\n');
fprintf(fid1,'-----  Model Configuration  ----------------------------------------------------\n');
fprintf(fid1,'%-21.0f',1);
fprintf(fid1,'NumSect:                   Number of circumferential sectors.\n');
fprintf(fid1,'%-21.0f',str2num(get(handles.MaxIter,'String')));
fprintf(fid1,'MaxIter:                   Max number of iterations for induction factor.\n');
fprintf(fid1,'%-21.1e',str2num(get(handles.ATol,'String')));
fprintf(fid1,'ATol:                      Error tolerance for induction iteration.\n');
fprintf(fid1,'%-21.1e',1.0e-006);
fprintf(fid1,'SWTol:                     Error tolerance for skewed-wake iteration.\n');
fprintf(fid1,'-----  Algorithm Configuration  ------------------------------------------------\n');
fprintf(fid1,'%-21s',handles.tiploss);
fprintf(fid1,'TipLoss:                   Use the Prandtl tip-loss model?\n');
fprintf(fid1,'%-21s',handles.hubloss);
fprintf(fid1,'HubLoss:                   Use the Prandtl hub-loss model?\n');
fprintf(fid1,'%-21s',handles.swirl);
fprintf(fid1,'Swirl:                     Include Swirl effects?\n');
fprintf(fid1,'%-21s','True ');
fprintf(fid1,'SkewWake:                  Apply skewed-wake correction?\n');
fprintf(fid1,'%-21s',handles.advbrake);
fprintf(fid1,'AdvBrake:                  Use the advanced brake-state model?\n');
fprintf(fid1,'%-21s',handles.indprop);
fprintf(fid1,'IndProp:                   Use PROP-PC instead of PROPX induction algorithm?\n');
fprintf(fid1,'%-21s',handles.aidrag);
fprintf(fid1,'AIDrag:                    Use the drag term in the axial induction calculation?\n');
fprintf(fid1,'%-21s',handles.tidrag);
fprintf(fid1,'TIDrag:                    Use the drag term in the tangential induction calculation?\n');
fprintf(fid1,'-----  Turbine Data  -----------------------------------------------------------\n');
fprintf(fid1,'%-21.0f',str2num(get(handles.NumBlade,'String')));
fprintf(fid1,'NumBlade:                  Number of blades.\n');
fprintf(fid1,'%-21.6f',0.5*str2num(get(handles.RotorDia,'String')));
fprintf(fid1,'RotorRad:                  Rotor radius [length].\n');
fprintf(fid1,'%-21.6f',0.5*str2num(get(handles.HubDia,'String')));
fprintf(fid1,'HubRad:                    Hub radius [length or div by radius].\n');
fprintf(fid1,'%-21.1f',str2num(get(handles.PreCone,'String')));
fprintf(fid1,'PreCone:                   Precone angle, positive downstream [deg].\n');
fprintf(fid1,'%-21.1f',str2num(get(handles.ShaftTilt,'String')));
fprintf(fid1,'Tilt:                      Shaft tilt [deg].\n');
fprintf(fid1,'%-21.1f',0);
fprintf(fid1,'Yaw:                       Yaw error [deg].\n');
fprintf(fid1,'%-21.6f',str2num(get(handles.HubHt,'String')));
fprintf(fid1,'HubHt:                     Hub height [length or div by radius].\n');
fprintf(fid1,'%-21.0f',str2num(get(handles.NumSeg,'String')));
fprintf(fid1,'NumSeg:                    Number of blade segments (entire rotor radius).\n');
fprintf(fid1,'RElm    Twist   Chord  AFfile  PrntElem\n');

%Build the 'Aerodynamic Inputs.wtp' text file
fprintf(fid2,'-----  Aerodynamic Data  -------------------------------------------------------\n');
fprintf(fid2,'%-21.6f',str2num(get(handles.Rho,'String')));
fprintf(fid2,'Rho:                 Air density [mass/volume].\n');
fprintf(fid2,'%-21.10f',str2num(get(handles.KinVisc,'String')));
fprintf(fid2,'KinVisc:             Kinematic air viscosity\n');
fprintf(fid2,'%-21.3f',0);
fprintf(fid2,'ShearExp:            Wind shear exponent (1/7 law = 0.143).\n');
    if(get(handles.CheckCavit,'Value'))
        fprintf(fid2,'%-21s','True ');
    else
        fprintf(fid2,'%-21s','False');
    end
fprintf(fid2,'UseCm:               Are Cm data included in the airfoil tables?\n');

%Build the "Output Configurations.wtp" text file
fprintf(fid3,'-----  I/O Settings  -----------------------------------------------------------\n');
fprintf(fid3,'%-21s','True ');
fprintf(fid3,'TabDel:                    Make output tab-delimited (fixed-width otherwise).\n');
fprintf(fid3,'%-21s','True ');
fprintf(fid3,'KFact:                     Output dimensional parameters in K (e.g., kN instead on N)\n');
fprintf(fid3,'%-21s','True ');
fprintf(fid3,'WriteBED:                  Write out blade element data to "<rootname>.bed"?\n');
fprintf(fid3,'%-21s','False');
fprintf(fid3,'InputTSR:                  Input speeds as TSRs?\n');
fprintf(fid3,'%-21s','"mps"');
fprintf(fid3,'SpdUnits:                  Wind-speed units (mps, fps, mph).\n');
fprintf(fid3,'-----  Combined-Case Analysis  -------------------------------------------------\n');
fprintf(fid3,'%-21.0f',0);
fprintf(fid3,'NumCases:                  Number of cases to run.  Enter zero for parametric analysis.\n');
fprintf(fid3,'WS or TSR   RotSpd   Pitch                      Remove following block of lines if NumCases is zero.\n');
fprintf(fid3,'-----  Parametric Analysis (Ignored if NumCases > 0 )  -------------------------\n');
fprintf(fid3,'%-21.0f',3);
fprintf(fid3,'ParRow:                    Row parameter    (1-rpm, 2-pitch, 3-tsr/speed).\n');
fprintf(fid3,'%-21.0f',1);
fprintf(fid3,'ParCol:                    Column parameter (1-rpm, 2-pitch, 3-tsr/speed).\n');
fprintf(fid3,'%-21.0f',2);
fprintf(fid3,'ParTab:                    Table parameter  (1-rpm, 2-pitch, 3-tsr/speed).\n');
fprintf(fid3,'%-21s','True ');
fprintf(fid3,'OutPwr:                    Request output of rotor power?\n');
fprintf(fid3,'%-21s','True ');
fprintf(fid3,'OutCp:                     Request output of Cp?\n');
fprintf(fid3,'%-21s','False');
fprintf(fid3,'OutTrq:                    Request output of shaft torque?\n');
fprintf(fid3,'%-21s','False');
fprintf(fid3,'OutFlp:                    Request output of flap bending moment?\n');
fprintf(fid3,'%-21s','False');
fprintf(fid3,'OutThr:                    Request output of rotor thrust?\n');
fprintf(fid3,'%-21s','0, 0, 0');
fprintf(fid3,'                           PitSt, PitEnd, PitDel:            First, last, delta blade pitch (deg).\n');
fprintf(fid3,'%-21s','0, 0, 0');
fprintf(fid3,'                           OmgSt, OmgEnd, OmgDel:            First, last, delta rotor speed (rpm).\n');
fprintf(fid3,'%-21s',[get(handles.SpdSt,'String') ', ' get(handles.SpdEnd,'String') ', ' get(handles.SpdDel,'String')]);
fprintf(fid3,'                          SpdSt, SpdEnd, SpdDel:            First, last, delta speeds.\n');

%close all open files
fclose('all');

guidata(hObject, handles);
%=========================================================================%

%=========================================================================%
%This section of the code executes when the "Begin Optimization" button is
%pressed.  This section performs the following tasks:
%          -reads the input variables from the GUI
%          -creates a new directory for the output files
%          -creates a new directory for the 3D Airfoil Data (if needed)
%          -reads in the initial population from "Initial_Population.xls"
%          -reads the custom flow probability from "Flow_Distribution.dat"
%           and interpolates the flow probability to the proper flow speeds
%          -defines all the Genetic Algorithm options
%          -checks the user input from the GUI for errors before
%           initializing the GA, then initializes the GA
%          -writes the best individual found by the GA back to the
%           "Initial_Population.xls" file

% --- Executes on button press in Optimize.
function Optimize_Callback(hObject, eventdata, handles)
% hObject    handle to Optimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
clear functions %reinitalizes the persistent variables

fprintf(1,'User inputs defined. Initializing...\n')

%Calls the "Save" function which writes all the WT_Perf template files
Save_Callback(hObject, eventdata, handles);

global SpdSt SpdEnd SpdDel OmgMin OmgMax Prated rho Patm Pv HubHt WatDepth Type...
       NumSeg HubRad Thickness_values ThickMethod RotorRad NumBlade family RootDir filename_main NumAFPoints...
       Correct_3D SpdCtrl PitCtrl CavSF OptMethod StructuralOpt ProbDist U_mean Weib_Umean Weib_k Weib_c user_pU_interp...
       ShaftTilt PreCone NumVars NumGens PopSize EliteCount ParetoFraction CrossFrc KinVisc CircleRoot...
       minRootChord maxRootChord RootTranSt RootTranEnd Normalized_AF_Coordinates...
       TwistLB TwistUB ChordLB ChordUB ThickLB ThickUB LB UB...
       Einput allowableStrain matDensity DecreasingST STmin STdel SFstruct...
       ParetoSet Pareto_filename OpenFile_Error RadialSpacing RecordFailures; %NumAFPoints added

%Read in the input variables from the GUI
SpdCtrl = str2num(handles.SpdCtrl);
PitCtrl = str2num(handles.PitCtrl);
NumBlade = str2num(get(handles.NumBlade,'String'));
NumSeg = str2num(get(handles.NumSeg,'String'));
RotorRad = 0.5*str2num(get(handles.RotorDia,'String'));
HubRad = 0.5*str2num(get(handles.HubDia,'String'));
HubHt = str2num(get(handles.HubHt,'String'));
Prated = str2num(get(handles.RatedPwr,'String'));
ShaftTilt = str2num(get(handles.ShaftTilt,'String'));
PreCone = str2num(get(handles.PreCone,'String'));
SpdSt = str2num(get(handles.SpdSt,'String'));
SpdEnd = str2num(get(handles.SpdEnd,'String'));
SpdDel = str2num(get(handles.SpdDel,'String'));
OmgMin = str2num(get(handles.OmgMin,'String'));
OmgMax = str2num(get(handles.OmgMax,'String'));
rho = str2num(get(handles.Rho,'String'));
KinVisc = str2num(get(handles.KinVisc,'String'));
CC = get(handles.CheckCavit,'Value');
Type = CC+1;
Pv = str2num(get(handles.VapPres,'String'));
Patm = str2num(get(handles.Patm,'String'));
WatDepth = str2num(get(handles.WatDepth,'String'));
CavSF = str2num(get(handles.CavSF,'String'));
Correct_3D = str2num(handles.dragmodel);
ThickMethod = str2num(handles.ThickMethod);
CircleRoot = get(handles.CircularRoot,'Value');
minRootChord = str2num(get(handles.minRtChord,'String'));
maxRootChord = str2num(get(handles.maxRtChord,'String'));
RootTranSt = str2num(get(handles.RtTranSt,'String'));
RootTranEnd = str2num(get(handles.RtTranEnd,'String'));

StructuralOpt = get(handles.StructuralOpt,'Value');
RecordFailures = get(handles.RecFail,'Value');

%blade element spacing
RadialSpacing = str2num(handles.Elm_spacing);

Einput = str2num(get(handles.Einput,'String')); %GPa bulk material elasticity
allowableStrain = str2num(get(handles.allowableStrain,'String'))./1e6; %user input in microstrain, convert to strain
matDensity = str2num(get(handles.matDensity,'String')); %(kg/m^3) bulk material density
STmin = str2num(get(handles.STmin,'String'))/1000; %minimum shell thickness (m)
STdel = str2num(get(handles.STdel,'String'))/1000; %shell thickness increment (m)
SFstruct = str2num(get(handles.SFstruct,'String')); %safety factor multiplied to bending moments
DecreasingST = str2num(handles.DecreasingST);

OptMethod = str2num(handles.OptMethod);
ProbDist = str2num(handles.ProbDist);
opened_file = get(handles.FlowDist_Filename,'String');    
    if strcmp(opened_file,'no file selected')
        OpenFile_Error = 1;
    else
        OpenFile_Error = 0;
    end

U_mean = str2num(get(handles.U_box,'String'));
Weib_Umean = str2num(get(handles.Weib_U_box,'String'));
Weib_k = str2num(get(handles.k_box,'String'));
Weib_c = str2num(get(handles.c_box,'String'));
    
NumGens = str2num(get(handles.NumGens,'String'));
PopSize = str2num(get(handles.PopSize,'String'));
EliteCount = str2num(get(handles.EliteCount,'String'));
ParetoFraction = str2num(get(handles.ParetoFraction,'String'));
CrossFrc = str2num(get(handles.CrossFrc,'String'));
GATol = str2num(get(handles.GATol,'String'));
TwistLB = [str2num(get(handles.TwLB1,'String')) str2num(get(handles.TwLB2,'String')) str2num(get(handles.TwLB3,'String')) str2num(get(handles.TwLB4,'String')) str2num(get(handles.TwLB5,'String'))];
TwistUB = [str2num(get(handles.TwUB1,'String')) str2num(get(handles.TwUB2,'String')) str2num(get(handles.TwUB3,'String')) str2num(get(handles.TwUB4,'String')) str2num(get(handles.TwUB5,'String'))];
ChordLB = [str2num(get(handles.CLB1,'String')) str2num(get(handles.CLB2,'String')) str2num(get(handles.CLB3,'String')) str2num(get(handles.CLB4,'String')) str2num(get(handles.CLB5,'String'))];
ChordUB = [str2num(get(handles.CUB1,'String')) str2num(get(handles.CUB2,'String')) str2num(get(handles.CUB3,'String')) str2num(get(handles.CUB4,'String')) str2num(get(handles.CUB5,'String'))];
ThickLB = str2num(get(handles.ThickLB,'String'));
ThickUB = str2num(get(handles.ThickUB,'String'));

%if user specified no bounds, or used any NaN values, set bounds automatically
if (isempty(ThickLB) || isempty(ThickUB)) && length(Thickness_values) > 1
    gaThickLB = zeros(1,length(Thickness_values));
    gaThickUB = ones(1,length(Thickness_values));
    if ThickMethod == 1
        gaThickLB = ThickLB(1:end-1);
        gaThickUB = ThickUB(1:end-1);
    end
elseif (isempty(ThickLB) || isempty(ThickUB)) && length(Thickness_values) == 1
    gaThickLB = [];
    gaThickUB = [];
elseif any(isnan(ThickLB)) || any(isnan(ThickUB))
    %replace the NaN values with appropriate values
    gaThickLB = ThickLB;
    gaThickUB = ThickUB;
    gaThickLB(isnan(gaThickLB)) = 0;
    gaThickUB(isnan(gaThickUB)) = 1;
elseif all(isfinite(ThickLB)) || all(isfinite(ThickUB))
    gaThickLB = ThickLB;
    gaThickUB = ThickUB;
end
LB = [TwistLB ChordLB gaThickLB]; %gets passed into ga function
UB = [TwistUB ChordUB gaThickUB]; %gets passed into ga function
         
family = get(handles.Family,'String');
Thickness_values = sort(str2num(get(handles.StepVals,'String')),'descend');  %Thickness values must be in descending order
filename_main = get(handles.Filename,'String');
    
%Create a new output directory for the user's filename, and copy the WT_Perf executable to this directory
RootDir = pwd; %pwd is a function which returns the present working directory

%Check if output directory already exists, if yes, then delete and recreate
if exist([RootDir '\Output_Files\' filename_main],'dir') == 7;
   rmdir([RootDir '\Output_Files\' filename_main],'s');
   mkdir([RootDir '\Output_Files\' filename_main]);
else
   mkdir([RootDir '\Output_Files\' filename_main]);
end

%save a screenshot of the GUI to output directory
saveas(gcf,[RootDir '\Output_Files\' filename_main '\' filename_main '_Input.bmp'],'bmp');

mkdir([RootDir '\Output_Files\' filename_main '\Airfoil_Data']);
copyfile([RootDir '\Input_Files\Templates\WT_Perf_Jan1.exe'],[RootDir '\Output_Files\' filename_main '\WT_Perf_Jan1.exe']);

%Copy the template for the Excel output file to the new output directory
copyfile([RootDir '\Input_Files\Templates\Output_Template.xls'],[RootDir '\Output_Files\' filename_main '\' filename_main  '_Output.xls']);

%If a Stall-Delay model is being used, need to create a directory for the 3D Airfoil Data as well
if Correct_3D ~= 0
mkdir([RootDir '\Output_Files\' filename_main '\Airfoil_Data\3D_Airfoil_Data'])
end

%If we are recording the failed cases, make new folder
if RecordFailures == 1
   mkdir([RootDir '\Output_Files\' filename_main '\Failed_Cases'])
   copyfile([RootDir '\Input_Files\Templates\WT_Perf_Jan1.exe'],[RootDir '\Output_Files\' filename_main '\Failed_Cases\WT_Perf_Jan1.exe']);
end



%Establish the Linear Inequality Constraints Ax<=b. Establish constraints
%such that the twist, chord, and %thickness distributions are monotonically decreasing
A = [-1 1 0 0 0 0 0 0 0 0;   %c1>=c2
      0 -1 1 0 0 0 0 0 0 0;  %c2>=c3
      0 0 -1 1 0 0 0 0 0 0;  %c3>=c4
      0 0 0 -1 1 0 0 0 0 0;  %c4>=c5
      0 0 0 0 0 -1 1 0 0 0;  %tw1>=tw2
      0 0 0 0 0 0 -1 1 0 0;  %tw2>=tw3
      0 0 0 0 0 0 0 -1 1 0;  %tw3>=tw4
      0 0 0 0 0 0 0 0 -1 1]; %tw4>=tw5

if length(Thickness_values) == 1; %only 1 airfoil profile exists
    NumVars = 10;
elseif length(Thickness_values) > 1 && ThickMethod == 1 %Piecewise Constant
    NumVars = 10 + length(Thickness_values)-1; %multiple airfoil profiles exist
    A(:,NumVars) = 0;
    nn=0;
    for n = 9:6+length(Thickness_values)
    A(n,11+nn) = 1;
    A(n,11+nn+1) = -1;
    nn=nn+1;
    end  
elseif length(Thickness_values) > 1 && ThickMethod == 2 %Piecewise Linear
    NumVars = 10 + length(Thickness_values); %multiple airfoil profiles exist
    nn=0;
    for n = 9:7+length(Thickness_values)
    A(n,11+nn) = 1;
    A(n,11+nn+1) = -1;
    nn=nn+1;
    end 
end 
b = zeros(size(A,1),1); %Linear Inequality Constraints Ax<=b

if SpdCtrl == 0   %Fixed Speed: rotor speed becomes the last variable
   NumVars = NumVars + 1; 
   LB = [LB OmgMin]; %add the min/max rotors speeds into the bounds
   UB = [UB OmgMax];
   A(:,NumVars) = 0;
end

    %=====================================================================%
%     NumSeedIndiv = xlsread('Input_Files\Initial_Population.xls',1,'B6');
%     SeedInitPop = get(handles.SeedInitPop,'Value');
%     
%     %Read in Initial Populations from the "Initial Population.xls" file
%     %NOTE: can only seed individuals with the same number of variables as
%     %the user specified in the current HARP_Opt session
%     
%     if SeedInitPop == 1 && NumSeedIndiv ~= 0;
%         InitPop = xlsread('Input_Files\Initial_Population.xls',1,['A10:P' num2str(10 + NumSeedIndiv - 1)]);
%         if NumVars == 16;       %Rotor Speed is the 16th variable       
%         EmptyCells = find(isnan(InitPop));
%         OmgAvg = (OmgMin + OmgMax)/2;
%         InitPop(EmptyCells) = OmgAvg;  
%             if size(InitPop,2) < 16; %the user has neglected to seed a 16th variable
%             InitPop(:,16) = OmgAvg;  %If the inital population does not contain a value for Rotor Spd., use this average value
%             end
%         else                       %Only need 15 variables
%         InitPop = InitPop(:,1:15); %get rid of the 16th variables
%         end
%     else
        InitPop = [];
%     end
    %=====================================================================%

%Configure the options for the Genetic Algorithm

%'PopInitRange',[LB;UB], ---
%{@gaplotbestindiv,@gaplotbestf,@gaplotdistance,
GAoptions = gaoptimset('InitialPopulation',InitPop,'CreationFcn',{@gaCustom_Creation},'Generations',NumGens,'PopulationSize',PopSize,'EliteCount',EliteCount,...
                     'CrossoverFraction',CrossFrc,'TolFun',GATol,'CrossoverFcn',{@crossoverintermediate},...
                     'SelectionFcn',{@selectionstochunif},'FitnessScalingFcn',{@fitscalingrank},'MutationFcn',{@mutationadaptfeasible},...
                     'PlotFcns',{@gaCustomPlot},'OutputFcns',@gaoutputfcn,'Display','diagnose');

MOGAoptions = gaoptimset('InitialPopulation',InitPop,'CreationFcn',{@gaCustom_Creation},'Generations',NumGens,'PopulationSize',PopSize,...
                     'CrossoverFraction',CrossFrc,'TolFun',GATol,'CrossoverFcn',{@crossoverintermediate},...
                     'MutationFcn',{@mutationadaptfeasible},'ParetoFraction',ParetoFraction,'PlotFcn',@gaCustomPlot,...
                     'OutputFcns',@gaoutputfcn,'Display','diagnose');
                 
%Need to check for user input errors before starting the Genetic Algorithm
UserInputError = Error_Check; %Calls the Error_Check.m function
if UserInputError == 0;
    
    %=====================================================================%
    %Read in variables and perform tasks related to AEP
       if ProbDist == 3;
           %User Defined Distribution: read in data from "Flow_Distribution.dat"
           fid = fopen(get(handles.Hidden_FlowDist_Filename,'String'),'rt');
           user_data = textscan(fid,'%f %f','HeaderLines',13);
           fclose(fid);
           U_custom = user_data{1};
           pU_custom = user_data{2};
           %Now interpolate the custom p(U) distribution
           V = (SpdSt:SpdDel:SpdEnd)';
           user_pU_interp = interp1(U_custom,pU_custom,V,'pchip');
           %if the spline caused any negative values change them to zero
           user_pU_interp(user_pU_interp<0) = 0;
       else
           user_pU_interp = [];    
       end  

    %=====================================================================%
 
    
    %If the airfoil files can be found (i.e. there were no user input errors,
    %then copy the 2D airfoil files that will be used to the new ouput directory
    if CircleRoot == 1;
        ThickVals = [100 Thickness_values];
    else
        ThickVals = Thickness_values;
    end    
    
    for n = 1:length(ThickVals)
        if ThickVals(n) >= 99.95;
        Thick_Suffix = '_1000';
        elseif ThickVals(n) >= 9.95 && ThickVals(n) < 99.95;
        Thick_Suffix = ['_0' num2str(10*ThickVals(n),'%3.0f')];
        else  
        Thick_Suffix = ['_00' num2str(10*ThickVals(n),'%3.0f')];
        end

%         AFfile(n,:) = [family Thick_Suffix '.dat'];
        AFfile(n,:) = [family Thick_Suffix];
        copyfile([RootDir '\Input_Files\Airfoil_Data\' AFfile(n,:) '.dat'],[RootDir '\Output_Files\' filename_main '\Airfoil_Data\' AFfile(n,:) '.dat']);
    end

    %We need to create new airfoil files which are interpolated between the available airfoil percent thicknesses
    thick_stepsize = 0.1; %airfoil thickness will be interpolated using this interval stepsize
    if CircleRoot == 0 && ThickMethod == 1
        %do nothing
    else
        for n = 1:length(ThickVals)-1;
        Interp_Thickness_Coefs([AFfile(n,:) '.dat'],[AFfile(n+1,:) '.dat'],ThickVals(n),ThickVals(n+1),thick_stepsize); 
        end    
    end
    
    if StructuralOpt == 1;
       mkdir([RootDir '\Output_Files\' filename_main '\Airfoil_Data\Coordinates']);
          
        if length(ThickVals) > 1

           for n = 1:length(ThickVals)-1;    
                [AFcoordinates] = Interp_Thickness_Profile([AFfile(n,:) '.prof'],[AFfile(n+1,:) '.prof'],ThickVals(n),ThickVals(n+1),thick_stepsize);
                if n == 1
                Normalized_AF_Coordinates = AFcoordinates;
                else
                Normalized_AF_Coordinates = [Normalized_AF_Coordinates;AFcoordinates(2:end,:)];
                end
           end
       elseif length(ThickVals) == 1;
            %only one airfoil being used
            in_path1 = fopen([RootDir '\Input_Files\Airfoil_Data\' AFfile '.prof'],'rt');
            A1 = cell2mat(textscan(in_path1,'%f %f','HeaderLines',4,'CollectOutput',1));
            fclose(in_path1);
            x1 = A1(:,1); %x coordinates
            y1 = A1(:,2); %y coordinates
            NumAFPoints = 40; %Number of points for new file, actual # of points will be NumAFPoints -1
            x_upper = hcosspace(0,1,NumAFPoints/2,3);
            x_lower = flipud(x_upper(1:(end-1)));
            Xnew = [x_upper;x_lower];
            a1 = find(A1(:,1)==1);
            y_upper1 = interp1(x1(1:a1),y1(1:a1),x_upper,'pchip');
            y_lower1 = interp1(x1(a1:end),y1(a1:end),x_lower,'pchip');
            Ynew1 = [y_upper1;y_lower1];
            Normalized_AF_Coordinates = {[Xnew Ynew1], ThickVals};
            %Now write the newly interpolated profile to the output directory
            out_path1 = fopen([RootDir '\Output_Files\' filename_main '\Airfoil_Data\Coordinates\' AFfile '.prof'],'Wt');
            %Print the headers and data for the newly interpolated files
            fprintf(out_path1,'This file was generated automatically by HARP_Opt.\n');
            fprintf(out_path1,'These airfoil coordinates were interpolated from the file %s to have %g points, using cosine spacing.\n',[AFfile '.prof'],NumAFPoints-1);
            fprintf(out_path1,'X\tY\n\n');
            fprintf(out_path1,'%3.6f\t%3.6f\n',[Xnew Ynew1]');
            fclose(out_path1);
        end
          coords_filename = [RootDir '\Output_Files\' filename_main '\Airfoil_Data\Coordinates\Normalized_AF_Coordinates'];
          save(coords_filename,'Normalized_AF_Coordinates')
    end
    
    if StructuralOpt == 1;
  
    [ParetoSet,Fvals,exitFlag,Output] = gamultiobj(@Main,NumVars,A,b,[],[],LB,UB,MOGAoptions);

    %sort the Pareto set
    [sortedFvals sortIndex] = sortrows(Fvals);
    sortedParetoSet = zeros(size(ParetoSet));
    for n = 1:length(sortIndex)
        sortedParetoSet(n,:) = ParetoSet((sortIndex(n)),:);
    end
    
    Fvals = sortedFvals;
    ParetoSet = sortedParetoSet;%pause on
    
    NaNcols = isnan(Fvals(:,1));
    Fvals(NaNcols,:) = [];
    ParetoSet(NaNcols,:) = [];
    
    %disp('ParetoSet');
    if size(ParetoSet,1) > 1
        for n = 1:size(ParetoSet,1)
             disp('ParetoSet');disp(n);
            disp(ParetoSet(n,:));
        end
    end
    %pause on;
    %keyboard;
    
    dumpfile = [RootDir '\Output_Files\' filename_main '\' filename_main  '_dump.txt'];
    save dumpfile cell2mat(ParetoSet) -ASCII
        
    %Start writing the Pareto data to Excel file
    Excel = actxserver('excel.application');
    %Excel.Visible = 1;
    ExcelFile = [RootDir '\Output_Files\' filename_main '\' filename_main  '_Output.xls'];
    % MS EXCEL 2007
    %ExcelFileTemplate = [RootDir '\Output_Files\' filename_main '\' filename_main  '_Output.xltx'];
    % MS EXCEL 2003
    %ExcelFileTemplate = [RootDir '\Output_Files\' filename_main '\' filename_main  '_Output.xlt'];
    %open the excel file, full path needs to be mentioned or else excel will pick it from most recently opened files.
    ExcelWorkbook = Excel.Workbooks.Open(ExcelFile); %This file has a sheet named "1" which has some formating
    wksheet = ExcelWorkbook.Worksheets.Item('1'); % Choose desired sheet
    
    
    
    %Make a copy of the Excel template sheet for the multiple Pareto solutions
    if size(ParetoSet,1) > 1
        for n = 2:size(ParetoSet,1)
            %make a copy of the template sheet % PROBLEM W/ EXCEL 2000,
            %2007.. does this fail on >19 sets?
            % copy failure recorded: see support microsoft.com/kb/210684
            wksheet.Copy(wksheet); %this will create a sheet called "1 (2)" and places it before "1" (not sure if this will be consistent)
            %wksheet.Add Type=ExcelFileTemplate
            newSheet=ExcelWorkbook.Worksheets.Item('1 (2)'); %get a handle to this copied sheet
            newSheet.Name=num2str(n); %rename it with a new name
            %if mod(n,10)==0
            %ExcelWorkbook.Save 
            %ExcelWorkbook.Close(false)  
            %ExcelWorkbook = Excel.Workbooks.Open(ExcelFile); 
            %wksheet = ExcelWorkbook.Worksheets.Item(num2str(n));
        end
        
        % Get the list of sheets in the workbook
        Sheets = ExcelWorkbook.Sheets;
        % Moves Sheet "1" to before sheet "2"
        Sheets.Item('1').Move(Sheets.Item('2'))          
    end
    
    ExcelWorkbook.Save % save the changes
    ExcelWorkbook.Close(false)  % Close Excel workbook.
    Excel.Quit % close excel
    delete(Excel); 
    %system('taskkill /F /IM EXCEL.EXE');
    
        %Ok, now we have copied all the sheets we needed, now let's go through
        %the sheets one by one and write the data
        for n = 1:size(ParetoSet,1)
            Pareto_filename = [filename_main '_' num2str(n)]; %this is a global variable which gets changed, Post_Process.m uses this variable
            Main([ParetoSet(n,:)';n]);
        end
    
    else
    %Initialze the Genetic Algorithm and return the best individual found
    [x_final] = ga(@Main,NumVars,A,b,[],[],LB,UB,[],GAoptions);
    
    Pareto_filename = filename_main; %this is a global variable which gets changed, Post_Process.m uses this variable
    ParetoSet = x_final;
    Main([x_final 1]);
    end
    


    %Append the best individual to the "Initial Population.xls" file
%     xlswrite('Input_Files\Initial_Population.xls',NumSeedIndiv+1,'B6:B6');
%     if NumVars == 15;
%     xlswrite('Input_Files\Initial_Population.xls',x_final',['A' num2str(10 + NumSeedIndiv) ':O' num2str(10 + NumSeedIndiv)]);
%     else
%     xlswrite('Input_Files\Initial_Population.xls',x_final',['A' num2str(10 + NumSeedIndiv) ':P' num2str(10 + NumSeedIndiv)]);
%     end

    disp('...Optimization Complete.'); %The entire program has finished running!
    if StructuralOpt == 1
    fprintf('The number of points on the Pareto front was: %d\n', size(Fvals,1));
    end
else
    disp('User input error. Program terminated.');
    fclose('all');
end
%=========================================================================%

%=========================================================================%
%The following sections of code are used to create objects in the GUI. Most
%of this code was automatically generated by the MATLAB GUI builder, but
%some modifications are made to this code to perform special actions, such
%as automatically calculating numbers & displaying them in the GUI, or for
%automatically "greying" out user inputs when appropriate.

% --- Executes during object creation, after setting all properties.
function GUI_Image_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GUI_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate GUI_Image

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Text Boxes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NumBlade_Callback(hObject, eventdata, handles)
% hObject    handle to NumBlade (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
 
%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function NumBlade_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumBlade (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NumSeg_Callback(hObject, eventdata, handles)
% hObject    handle to HubHt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
 
%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function NumSeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HubHt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RotorDia_Callback(hObject, eventdata, handles)
% hObject    handle to NumBlade (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
 
RPM = str2num(get(handles.OmgMax,'String'));
R = 0.5*str2num(get(handles.RotorDia,'String'));
set(handles.TipSpd,'String',num2str(RPM*R*pi/30,'%6.1f'));

hr = str2num(get(handles.HubDia,'String'));
rr = str2num(get(handles.RotorDia,'String'));
if str2num(get(handles.RtTranSt,'String')) < hr/rr;
set(handles.RtTranSt,'String',num2str(hr/rr))
end

%Define the radial locations of the control points: using hafl-cosine_spacing spacing
   if get(handles.CircularRoot,'Value') == 0
   rad_CP = hcosspace(0.5*str2num(get(handles.HubDia,'String')),0.5*str2num(get(handles.RotorDia,'String')),5,0);    
   else
   rad_CP = hcosspace(str2num(get(handles.RtTranEnd,'String'))*0.5*str2num(get(handles.RotorDia,'String')),0.5*str2num(get(handles.RotorDia,'String')),5,0);      
   end
set(handles.radCP1,'String',num2str(rad_CP(1),4));
set(handles.radCP2,'String',num2str(rad_CP(2),4));
set(handles.radCP3,'String',num2str(rad_CP(3),4));
set(handles.radCP4,'String',num2str(rad_CP(4),4));
set(handles.radCP5,'String',num2str(rad_CP(5),4));

UpdateEstimates(handles) %calls the function to estimate AEP and CF

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function RotorDia_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumBlade (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function HubDia_Callback(hObject, eventdata, handles)
% hObject    handle to HubDia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));

hr = str2num(get(handles.HubDia,'String'));
rr = str2num(get(handles.RotorDia,'String'));
if str2num(get(handles.RtTranSt,'String')) < hr/rr;
set(handles.RtTranSt,'String',num2str(hr/rr))
end

%Define the radial locations of the control points: using hafl-cosine_spacing spacing
   if get(handles.CircularRoot,'Value') == 0
   rad_CP = hcosspace(0.5*str2num(get(handles.HubDia,'String')),0.5*str2num(get(handles.RotorDia,'String')),5,0);    
   else
   rad_CP = hcosspace(str2num(get(handles.RtTranEnd,'String'))*0.5*str2num(get(handles.RotorDia,'String')),0.5*str2num(get(handles.RotorDia,'String')),5,0);      
   end
set(handles.radCP1,'String',num2str(rad_CP(1),4));
set(handles.radCP2,'String',num2str(rad_CP(2),4));
set(handles.radCP3,'String',num2str(rad_CP(3),4));
set(handles.radCP4,'String',num2str(rad_CP(4),4));
set(handles.radCP5,'String',num2str(rad_CP(5),4));
 
%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function HubDia_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HubDia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function HubHt_Callback(hObject, eventdata, handles)
% hObject    handle to HubHt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
 
%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function HubHt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HubHt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RatedPwr_Callback(hObject, eventdata, handles)
% hObject    handle to RatedPwr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
 
UpdateEstimates(handles) %calls the function to estimate AEP and CF

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function RatedPwr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RatedPwr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SpdSt_Callback(hObject, eventdata, handles)
% hObject    handle to SpdSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
 
UpdateEstimates(handles) %calls the function to estimate AEP and CF

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SpdSt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpdSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SpdEnd_Callback(hObject, eventdata, handles)
% hObject    handle to SpdEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
 
UpdateEstimates(handles) %calls the function to estimate AEP and CF

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SpdEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpdEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SpdDel_Callback(hObject, eventdata, handles)
% hObject    handle to SpdDel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));

UpdateEstimates(handles) %calls the function to estimate AEP and CF

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SpdDel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpdDel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function OmgMin_Callback(hObject, eventdata, handles)
% hObject    handle to OmgMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OmgMin as text
%        str2double(get(hObject,'String')) returns contents of OmgMin as a double
%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
 
%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function OmgMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OmgMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function OmgMax_Callback(hObject, eventdata, handles)
% hObject    handle to OmgMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OmgMax as text
%        str2double(get(hObject,'String')) returns contents of OmgMax as a double
%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));

RPM = str2num(get(handles.OmgMax,'String'));
R = 0.5*str2num(get(handles.RotorDia,'String'));
set(handles.TipSpd,'String',num2str(RPM*R*pi/30,'%6.1f'));
 
%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function OmgMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OmgMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TipSpd_Callback(hObject, eventdata, handles)
% hObject    handle to TipSpd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TipSpd as text
%        str2double(get(hObject,'String')) returns contents of TipSpd as a double
%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function TipSpd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TipSpd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Rho_Callback(hObject, eventdata, handles)
% hObject    handle to Rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));

UpdateEstimates(handles) %calls the function to estimate AEP and CF

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Rho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function KinVisc_Callback(hObject, eventdata, handles)
% hObject    handle to KinVisc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
 
%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function KinVisc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KinVisc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function VapPres_Callback(hObject, eventdata, handles)
% hObject    handle to VapPres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VapPres as text
%        str2double(get(hObject,'String')) returns contents of VapPres as a double

% --- Executes during object creation, after setting all properties.
function VapPres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VapPres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Patm_Callback(hObject, eventdata, handles)
% hObject    handle to Patm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Patm as text
%        str2double(get(hObject,'String')) returns contents of Patm as a double

% --- Executes during object creation, after setting all properties.
function Patm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Patm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function WatDepth_Callback(hObject, eventdata, handles)
% hObject    handle to WatDepth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WatDepth as text
%        str2double(get(hObject,'String')) returns contents of WatDepth as a double

% --- Executes during object creation, after setting all properties.
function WatDepth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WatDepth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CavSF_Callback(hObject, eventdata, handles)
% hObject    handle to CavSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CavSF as text
%        str2double(get(hObject,'String')) returns contents of CavSF as a double

% --- Executes during object creation, after setting all properties.
function CavSF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CavSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function U_box_Callback(hObject, eventdata, handles)
% hObject    handle to U_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UpdateEstimates(handles) %calls the function to estimate AEP and CF

% Hints: get(hObject,'String') returns contents of U_box as text
%        str2double(get(hObject,'String')) returns contents of U_box as a double

% --- Executes during object creation, after setting all properties.
function U_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to U_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function k_box_Callback(hObject, eventdata, handles)
% hObject    handle to k_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k_box as text
%        str2double(get(hObject,'String')) returns contents of k_box as a double
cc = str2num(get(handles.c_box,'String'));
kk = str2num(get(handles.k_box,'String'));
Weib_Umean = cc*gamma(1+1/kk);
set(handles.Weib_U_box,'String',num2str(Weib_Umean,'%6.2f'));

UpdateEstimates(handles) %calls the function to estimate AEP and CF

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function k_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function c_box_Callback(hObject, eventdata, handles)
% hObject    handle to c_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c_box as text
%        str2double(get(hObject,'String')) returns contents of c_box as a double
cc = str2num(get(handles.c_box,'String'));
kk = str2num(get(handles.k_box,'String'));
Weib_Umean = cc*gamma(1+1/kk);
set(handles.Weib_U_box,'String',num2str(Weib_Umean,'%6.2f'));

UpdateEstimates(handles) %calls the function to estimate AEP and CF

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function c_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Weib_U_box_Callback(hObject, eventdata, handles)
% hObject    handle to Weib_U_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
 
%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Weib_U_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Weib_U_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MaxIter_Callback(hObject, eventdata, handles)
% hObject    handle to MaxIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
 
%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MaxIter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ATol_Callback(hObject, eventdata, handles)
% hObject    handle to ATol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%store the contents of input1_editText as a string. if the string 
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
 
%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ATol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ATol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NumGens_Callback(hObject, eventdata, handles)
% hObject    handle to NumGens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumGens as text
%        str2double(get(hObject,'String')) returns contents of NumGens as a double

% --- Executes during object creation, after setting all properties.
function NumGens_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumGens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PopSize_Callback(hObject, eventdata, handles)
% hObject    handle to PopSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PopSize as text
%        str2double(get(hObject,'String')) returns contents of PopSize as a double

% --- Executes during object creation, after setting all properties.
function PopSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function EliteCount_Callback(hObject, eventdata, handles)
% hObject    handle to EliteCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EliteCount as text
%        str2double(get(hObject,'String')) returns contents of EliteCount as a double

% --- Executes during object creation, after setting all properties.
function EliteCount_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EliteCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CrossFrc_Callback(hObject, eventdata, handles)
% hObject    handle to CrossFrc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CrossFrc as text
%        str2double(get(hObject,'String')) returns contents of CrossFrc as a double

% --- Executes during object creation, after setting all properties.
function CrossFrc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CrossFrc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function GATol_Callback(hObject, eventdata, handles)
% hObject    handle to GATol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GATol as text
%        str2double(get(hObject,'String')) returns contents of GATol as a double

% --- Executes during object creation, after setting all properties.
function GATol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GATol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Family_Callback(hObject, eventdata, handles)
% hObject    handle to Family (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Family as text
%        str2double(get(hObject,'String')) returns contents of Family as a double

% --- Executes during object creation, after setting all properties.
function Family_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Family (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Filename_Callback(hObject, eventdata, handles)
% hObject    handle to Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Filename as text
%        str2double(get(hObject,'String')) returns contents of Filename as a double

% --- Executes during object creation, after setting all properties.
function Filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Radio Buttons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TipLoss_SelectionChangeFcn(hObject, eventdata)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.tiploss=get(eventdata.NewValue,'String');
guidata(hObject,handles)

function HubLoss_SelectionChangeFcn(hObject, eventdata)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.hubloss=get(eventdata.NewValue,'String');
guidata(hObject,handles)

function Swirl_SelectionChangeFcn(hObject, eventdata)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.swirl=get(eventdata.NewValue,'String');
guidata(hObject,handles)

function AdvBrake_SelectionChangeFcn(hObject, eventdata)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.advbrake=get(eventdata.NewValue,'String');
guidata(hObject,handles)

function IndProp_SelectionChangeFcn(hObject, eventdata)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.indprop=get(eventdata.NewValue,'String');
guidata(hObject,handles)

function AIDrag_SelectionChangeFcn(hObject, eventdata)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.aidrag=get(eventdata.NewValue,'String');
guidata(hObject,handles)

function TIDrag_SelectionChangeFcn(hObject, eventdata)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.tidrag=get(eventdata.NewValue,'String');
guidata(hObject,handles)

function DecST_SelectionChangeFcn(hObject, eventdata)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'True12'
        handles.DecreasingST = '1';
    case 'False12'
        handles.DecreasingST = '0';
    otherwise       
end
%updates the handles structure
guidata(hObject,handles)

%Thickness Distribution Button Group
set(handles.ThickDist,'SelectionChangeFcn',@ThickDist_SelectionChangeFcn);
guidata(hObject, handles);

function ThickDist_SelectionChangeFcn(hObject, eventdata)
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'ThickMethod1'
      %execute this code when ThickMethod1 is selected
      handles.ThickMethod='1';
 
    case 'ThickMethod2'
      %execute this code when ThickMethod2 is selected
      handles.ThickMethod='2';
 
    otherwise
       % Code for when there is no match.
 
end
%updates the handles structure
guidata(hObject, handles);

    %=====================================================================%
    % Stall Delay Model Button Group
function Correct_3D_SelectionChangeFcn(hObject, eventdata)
 
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'Correct_3D_0'
      %execute this code when Correct_3D_0 is selected
      handles.dragmodel='0';
 
    case 'Correct_3D_1'
      %execute this code when Correct_3D_1 is selected
      handles.dragmodel='1';
 
    case 'Correct_3D_2'
      %execute this code when Correct_3D_2 is selected  
      handles.dragmodel='2';
    otherwise
       % Code for when there is no match.
 
end
%updates the handles structure
guidata(hObject, handles);
    %=====================================================================%

    %=====================================================================%
    % Speed Control Button Group
function RotSpdCtrl_SelectionChangeFcn(hObject, eventdata)
 
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'SpdCtrl0'
      %execute this code when SpdCtrl0 is selected
      handles.SpdCtrl='0';
      set(handles.Correct_3D_1,'Enable','on');
      set(handles.Correct_3D_2,'Enable','on');
 
    case 'SpdCtrl1'
      %execute this code when SpdCtrl1 is selected
      handles.SpdCtrl='1';
      set(handles.Correct_3D, 'SelectedObject', handles.Correct_3D_0);
      handles.dragmodel='0';
      set(handles.Correct_3D_1,'Enable','off');
      set(handles.Correct_3D_2,'Enable','off');
    
    otherwise
       % Code for when there is no match.
 
end
%updates the handles structure
guidata(hObject, handles);
    %=====================================================================%

    %=====================================================================%
    % Pitch Control Button Group
function BldPitCtrl_SelectionChangeFcn(hObject, eventdata)
 
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'PitCtrl0'
      %execute this code when PitCtrl0 is selected
      handles.PitCtrl='0';
    case 'PitCtrl1'
      %execute this code when PitCtrl1 is selected
      handles.PitCtrl='1';
    case 'PitCtrl2'
      %execute this code when PitCtrl2 is selected
      handles.PitCtrl='2';
    otherwise
       % Code for when there is no match.
 
end
%updates the handles structure
guidata(hObject, handles);
    %=====================================================================%

    %=====================================================================%
    % --- Executes on button press in CheckCavit.
function CheckCavit_Callback(hObject, eventdata, handles)
% hObject    handle to CheckCavit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of CheckCavit
    if (get(hObject,'Value') == get(hObject,'Max'))
       % Checkbox is checked-take approriate action
       set(handles.VapPres,'Enable','on');
       set(handles.WatDepth,'Enable','on');
       set(handles.Patm,'Enable','on');
       set(handles.CavSF,'Enable','on');
    else
       % Checkbox is not checked-take approriate action
       set(handles.VapPres,'Enable','off');
       set(handles.WatDepth,'Enable','off');
       set(handles.Patm,'Enable','off');
       set(handles.CavSF,'Enable','off');
    end

guidata(hObject, handles);
    %=====================================================================%

    
% --- Executes on button press in StructuralOpt.
function StructuralOpt_Callback(hObject, eventdata, handles)
% hObject    handle to StructuralOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of StructuralOpt
    if (get(hObject,'Value') == get(hObject,'Max'))
       % Checkbox is checked-take approriate action
       set(handles.Einput,'Enable','on');
       set(handles.allowableStrain,'Enable','on');
       set(handles.matDensity,'Enable','on');
       set(handles.SFstruct,'Enable','on');
       set(handles.STmin,'Enable','on');
       set(handles.STdel,'Enable','on');
       set(handles.True12,'Enable','on');
       set(handles.False12,'Enable','on');
       set(handles.EliteCount,'Enable','off');
       set(handles.ParetoFraction,'Enable','on');
    else
       % Checkbox is not checked-take approriate action
       set(handles.Einput,'Enable','off');
       set(handles.allowableStrain,'Enable','off');
       set(handles.matDensity,'Enable','off');
       set(handles.SFstruct,'Enable','off');
       set(handles.STmin,'Enable','off');
       set(handles.STdel,'Enable','off');
       set(handles.True12,'Enable','off');
       set(handles.False12,'Enable','off');
       set(handles.EliteCount,'Enable','on');
       set(handles.ParetoFraction,'Enable','off');
    end

guidata(hObject, handles);


% --- Executes on button press in CircularRoot.
function CircularRoot_Callback(hObject, eventdata, handles)
% hObject    handle to CircularRoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if (get(hObject,'Value') == get(hObject,'Max'))
        set(handles.minRtChord,'Enable','on');
        set(handles.maxRtChord,'Enable','on');
        set(handles.RtTranSt,'Enable','on');
        set(handles.RtTranEnd,'Enable','on');
    else
        set(handles.minRtChord,'Enable','off');
        set(handles.maxRtChord,'Enable','off');
        set(handles.RtTranSt,'Enable','off');
        set(handles.RtTranEnd,'Enable','off');
    end
    
hr = str2num(get(handles.HubDia,'String'));
rr = str2num(get(handles.RotorDia,'String'));
if str2num(get(handles.RtTranSt,'String')) < hr/rr;
set(handles.RtTranSt,'String',num2str(hr/rr))
end

   %Define the radial locations of the control points: using hafl-cosine_spacing spacing
   if get(handles.CircularRoot,'Value') == 0
   rad_CP = hcosspace(0.5*str2num(get(handles.HubDia,'String')),0.5*str2num(get(handles.RotorDia,'String')),5,0);    
   else
   rad_CP = hcosspace(str2num(get(handles.RtTranEnd,'String'))*0.5*str2num(get(handles.RotorDia,'String')),0.5*str2num(get(handles.RotorDia,'String')),5,0);      
   end
set(handles.radCP1,'String',num2str(rad_CP(1),4));
set(handles.radCP2,'String',num2str(rad_CP(2),4));
set(handles.radCP3,'String',num2str(rad_CP(3),4));
set(handles.radCP4,'String',num2str(rad_CP(4),4));
set(handles.radCP5,'String',num2str(rad_CP(5),4));

% Hint: get(hObject,'Value') returns toggle state of CircularRoot
guidata(hObject, handles);

    %=====================================================================%
    % Flow Distribution Button Group
function ProbType_SelectionChangeFcn(hObject, eventdata)
 
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'ProbDist1'
      %execute this code when ProbDist1 (Rayleigh) is selected
      handles.ProbDist='1';
      set(handles.U_box,'Enable','on');
      set(handles.k_box,'Enable','off');
      set(handles.c_box,'Enable','off');
      set(handles.Weib_U_box,'Enable','off');
      set(handles.FlowFile_Button,'Enable','off');
      UpdateEstimates(handles) %calls the function to estimate AEP and CF

    case 'ProbDist2'
      %execute this code when ProbDist2 (Weibull) is selected
      handles.ProbDist='2';
      set(handles.U_box,'Enable','off');
      set(handles.k_box,'Enable','on');
      set(handles.c_box,'Enable','on');
      set(handles.Weib_U_box,'Enable','inactive');      
      set(handles.FlowFile_Button,'Enable','off');
      UpdateEstimates(handles) %calls the function to estimate AEP and CF

    case 'ProbDist3'
      %execute this code when ProbDist3 (custom) is selected  
      handles.ProbDist='3';
      set(handles.U_box,'Enable','off');
      set(handles.k_box,'Enable','off');
      set(handles.c_box,'Enable','off');
      set(handles.Weib_U_box,'Enable','off');
      set(handles.FlowFile_Button,'Enable','on');
      UpdateEstimates(handles) %calls the function to estimate AEP and CF

    case 'ProbDist4'
      %execute this code when ProbDist4 (none) is selected  
      handles.ProbDist='0';
      set(handles.U_box,'Enable','off');
      set(handles.k_box,'Enable','off');
      set(handles.c_box,'Enable','off');
      set(handles.Weib_U_box,'Enable','off');
      set(handles.FlowFile_Button,'Enable','off');
      UpdateEstimates(handles) %calls the function to estimate AEP and CF

    otherwise
      % Code for when there is no match.
      handles.ProbDist='0';
 
end
%updates the handles structure
guidata(hObject, handles);
    %=====================================================================%

    %=====================================================================%
    % Optimization Method Button Group
function OptimMethod_SelectionChangeFcn(hObject, eventdata)

%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'OptAEP'
      %execute this code when OptAEP is selected
      handles.OptMethod='1';
        set(handles.ProbDist1,'Enable','on');
        set(handles.ProbDist2,'Enable','on');
        set(handles.ProbDist3,'Enable','on');
        set(handles.FlowFile_Button,'Enable','off');        
        set(handles.ProbDist4,'Enable','off');
        set(handles.ProbType, 'SelectedObject', handles.ProbDist1);
        handles.ProbDist='1';
        set(handles.U_box,'Enable','on');
        set(handles.k_box,'Enable','off');
        set(handles.c_box,'Enable','off');
        set(handles.Weib_U_box,'Enable','off');
        UpdateEstimates(handles) %calls the function to estimate AEP and CF

    case 'OptEff'
      %execute this code when OptEff is selected
      handles.OptMethod='0';
      set(handles.ProbDist1,'Enable','on');
      set(handles.ProbDist2,'Enable','on');
      set(handles.ProbDist3,'Enable','on');
      set(handles.FlowFile_Button,'Enable','off');      
      set(handles.ProbDist4,'Enable','on');
      set(handles.U_box,'Enable','off');
      set(handles.k_box,'Enable','off');
      set(handles.c_box,'Enable','off');
      set(handles.Weib_U_box,'Enable','off');
      set(handles.ProbType, 'SelectedObject', handles.ProbDist4);
      handles.ProbDist='0';
      UpdateEstimates(handles) %calls the function to estimate AEP and CF

    otherwise
       % Code for when there is no match.
 
end
%updates the handles structure
guidata(hObject, handles);
    %=====================================================================%
    
    %=====================================================================%
% --- Executes on button press in SeedInitPop.
function SeedInitPop_Callback(hObject, eventdata, handles)
% hObject    handle to SeedInitPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SeedInitPop
    %=====================================================================%

    %=====================================================================%
    %Bounds
function TwLB1_Callback(hObject, eventdata, handles)
% hObject    handle to TwLB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TwLB1 as text
%        str2double(get(hObject,'String')) returns contents of TwLB1 as a double

% --- Executes during object creation, after setting all properties.
function TwLB1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TwLB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TwLB2_Callback(hObject, eventdata, handles)
% hObject    handle to TwLB2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TwLB2 as text
%        str2double(get(hObject,'String')) returns contents of TwLB2 as a double

% --- Executes during object creation, after setting all properties.
function TwLB2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TwLB2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TwLB3_Callback(hObject, eventdata, handles)
% hObject    handle to TwLB3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TwLB3 as text
%        str2double(get(hObject,'String')) returns contents of TwLB3 as a double

% --- Executes during object creation, after setting all properties.
function TwLB3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TwLB3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TwLB4_Callback(hObject, eventdata, handles)
% hObject    handle to TwLB4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TwLB4 as text
%        str2double(get(hObject,'String')) returns contents of TwLB4 as a double

% --- Executes during object creation, after setting all properties.
function TwLB4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TwLB4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TwLB5_Callback(hObject, eventdata, handles)
% hObject    handle to TwLB5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TwLB5 as text
%        str2double(get(hObject,'String')) returns contents of TwLB5 as a double

% --- Executes during object creation, after setting all properties.
function TwLB5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TwLB5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TwUB1_Callback(hObject, eventdata, handles)
% hObject    handle to TwUB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TwUB1 as text
%        str2double(get(hObject,'String')) returns contents of TwUB1 as a double

% --- Executes during object creation, after setting all properties.
function TwUB1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TwUB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TwUB2_Callback(hObject, eventdata, handles)
% hObject    handle to TwUB2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TwUB2 as text
%        str2double(get(hObject,'String')) returns contents of TwUB2 as a double

% --- Executes during object creation, after setting all properties.
function TwUB2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TwUB2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TwUB3_Callback(hObject, eventdata, handles)
% hObject    handle to TwUB3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TwUB3 as text
%        str2double(get(hObject,'String')) returns contents of TwUB3 as a double

% --- Executes during object creation, after setting all properties.
function TwUB3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TwUB3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TwUB4_Callback(hObject, eventdata, handles)
% hObject    handle to TwUB4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TwUB4 as text
%        str2double(get(hObject,'String')) returns contents of TwUB4 as a double

% --- Executes during object creation, after setting all properties.
function TwUB4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TwUB4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TwUB5_Callback(hObject, eventdata, handles)
% hObject    handle to TwUB5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TwUB5 as text
%        str2double(get(hObject,'String')) returns contents of TwUB5 as a double

% --- Executes during object creation, after setting all properties.
function TwUB5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TwUB5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CLB1_Callback(hObject, eventdata, handles)
% hObject    handle to CLB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CLB1 as text
%        str2double(get(hObject,'String')) returns contents of CLB1 as a double

% --- Executes during object creation, after setting all properties.
function CLB1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CLB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CLB2_Callback(hObject, eventdata, handles)
% hObject    handle to CLB2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CLB2 as text
%        str2double(get(hObject,'String')) returns contents of CLB2 as a double

% --- Executes during object creation, after setting all properties.
function CLB2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CLB2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CLB3_Callback(hObject, eventdata, handles)
% hObject    handle to CLB3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CLB3 as text
%        str2double(get(hObject,'String')) returns contents of CLB3 as a double

% --- Executes during object creation, after setting all properties.
function CLB3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CLB3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CLB4_Callback(hObject, eventdata, handles)
% hObject    handle to CLB4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CLB4 as text
%        str2double(get(hObject,'String')) returns contents of CLB4 as a double

% --- Executes during object creation, after setting all properties.
function CLB4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CLB4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CLB5_Callback(hObject, eventdata, handles)
% hObject    handle to CLB5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CLB5 as text
%        str2double(get(hObject,'String')) returns contents of CLB5 as a double

% --- Executes during object creation, after setting all properties.
function CLB5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CLB5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CUB1_Callback(hObject, eventdata, handles)
% hObject    handle to CUB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CUB1 as text
%        str2double(get(hObject,'String')) returns contents of CUB1 as a double

% --- Executes during object creation, after setting all properties.
function CUB1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CUB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CUB2_Callback(hObject, eventdata, handles)
% hObject    handle to CUB2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CUB2 as text
%        str2double(get(hObject,'String')) returns contents of CUB2 as a double

% --- Executes during object creation, after setting all properties.
function CUB2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CUB2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CUB3_Callback(hObject, eventdata, handles)
% hObject    handle to CUB3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CUB3 as text
%        str2double(get(hObject,'String')) returns contents of CUB3 as a double

% --- Executes during object creation, after setting all properties.
function CUB3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CUB3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CUB4_Callback(hObject, eventdata, handles)
% hObject    handle to CUB4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CUB4 as text
%        str2double(get(hObject,'String')) returns contents of CUB4 as a double

% --- Executes during object creation, after setting all properties.
function CUB4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CUB4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CUB5_Callback(hObject, eventdata, handles)
% hObject    handle to CUB5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CUB5 as text
%        str2double(get(hObject,'String')) returns contents of CUB5 as a double

% --- Executes during object creation, after setting all properties.
function CUB5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CUB5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function StepVals_Callback(hObject, eventdata, handles)
% hObject    handle to StepVals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StepVals as text
%        str2double(get(hObject,'String')) returns contents of StepVals as a double

% --- Executes during object creation, after setting all properties.
function StepVals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StepVals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ThickLB_Callback(hObject, eventdata, handles)
% hObject    handle to ThickLB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThickLB as text
%        str2double(get(hObject,'String')) returns contents of ThickLB as a double

% --- Executes during object creation, after setting all properties.
function ThickLB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThickLB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ThickUB_Callback(hObject, eventdata, handles)
% hObject    handle to ThickUB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThickUB as text
%        str2double(get(hObject,'String')) returns contents of ThickUB as a double

% --- Executes during object creation, after setting all properties.
function ThickUB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThickUB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function minRtChord_Callback(hObject, eventdata, handles)
% hObject    handle to minRtChord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minRtChord as text
%        str2double(get(hObject,'String')) returns contents of minRtChord as a double


% --- Executes during object creation, after setting all properties.
function minRtChord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minRtChord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxRtChord_Callback(hObject, eventdata, handles)
% hObject    handle to maxRtChord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxRtChord as text
%        str2double(get(hObject,'String')) returns contents of maxRtChord as a double

% --- Executes during object creation, after setting all properties.
function maxRtChord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxRtChord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RtTranSt_Callback(hObject, eventdata, handles)
% hObject    handle to RtTranSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hr = str2num(get(handles.HubDia,'String'));
rr = str2num(get(handles.RotorDia,'String'));
if str2num(get(handles.RtTranSt,'String')) < hr/rr;
set(handles.RtTranSt,'String',num2str(hr/rr))
end

if str2num(get(handles.RtTranEnd,'String')) <= str2num(get(handles.RtTranSt,'String'));
set(handles.RtTranEnd,'String',num2str(str2num(get(handles.RtTranSt,'String'))))
end

%Define the radial locations of the control points: using hafl-cosine_spacing spacing
   if get(handles.CircularRoot,'Value') == 0
   rad_CP = hcosspace(0.5*str2num(get(handles.HubDia,'String')),0.5*str2num(get(handles.RotorDia,'String')),5,0);    
   else
   rad_CP = hcosspace(str2num(get(handles.RtTranEnd,'String'))*0.5*str2num(get(handles.RotorDia,'String')),0.5*str2num(get(handles.RotorDia,'String')),5,0);      
   end
set(handles.radCP1,'String',num2str(rad_CP(1),4));
set(handles.radCP2,'String',num2str(rad_CP(2),4));
set(handles.radCP3,'String',num2str(rad_CP(3),4));
set(handles.radCP4,'String',num2str(rad_CP(4),4));
set(handles.radCP5,'String',num2str(rad_CP(5),4));

% Hints: get(hObject,'String') returns contents of RtTranSt as text
%        str2double(get(hObject,'String')) returns contents of RtTranSt as a double
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function RtTranSt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RtTranSt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RtTranEnd_Callback(hObject, eventdata, handles)
% hObject    handle to RtTranEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if str2num(get(handles.RtTranEnd,'String')) <= str2num(get(handles.RtTranSt,'String'));
set(handles.RtTranEnd,'String',num2str(str2num(get(handles.RtTranSt,'String'))))
end

%Define the radial locations of the control points: using hafl-cosine_spacing spacing
   if get(handles.CircularRoot,'Value') == 0
   rad_CP = hcosspace(0.5*str2num(get(handles.HubDia,'String')),0.5*str2num(get(handles.RotorDia,'String')),5,0);    
   else
   rad_CP = hcosspace(str2num(get(handles.RtTranEnd,'String'))*0.5*str2num(get(handles.RotorDia,'String')),0.5*str2num(get(handles.RotorDia,'String')),5,0);      
   end
set(handles.radCP1,'String',num2str(rad_CP(1),4));
set(handles.radCP2,'String',num2str(rad_CP(2),4));
set(handles.radCP3,'String',num2str(rad_CP(3),4));
set(handles.radCP4,'String',num2str(rad_CP(4),4));
set(handles.radCP5,'String',num2str(rad_CP(5),4));

% Hints: get(hObject,'String') returns contents of RtTranEnd as text
%        str2double(get(hObject,'String')) returns contents of RtTranEnd as a double
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function RtTranEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RtTranEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function radCP1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radCP1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function radCP2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radCP2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function radCP3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radCP3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function radCP4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radCP4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function radCP5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radCP5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function FlowDist_Filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radCP1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function PreCone_Callback(hObject, eventdata, handles)
% hObject    handle to PreCone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PreCone as text
%        str2double(get(hObject,'String')) returns contents of PreCone as a double

% --- Executes during object creation, after setting all properties.
function PreCone_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PreCone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ShaftTilt_Callback(hObject, eventdata, handles)
% hObject    handle to ShaftTilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ShaftTilt as text
%        str2double(get(hObject,'String')) returns contents of ShaftTilt as a double

% --- Executes during object creation, after setting all properties.
function ShaftTilt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ShaftTilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Einput_Callback(hObject, eventdata, handles)
% hObject    handle to Einput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Einput as text
%        str2double(get(hObject,'String')) returns contents of Einput as a double

% --- Executes during object creation, after setting all properties.
function Einput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Einput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function allowableStrain_Callback(hObject, eventdata, handles)
% hObject    handle to allowableStrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of allowableStrain as text
%        str2double(get(hObject,'String')) returns contents of allowableStrain as a double

% --- Executes during object creation, after setting all properties.
function allowableStrain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowableStrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function matDensity_Callback(hObject, eventdata, handles)
% hObject    handle to matDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of matDensity as text
%        str2double(get(hObject,'String')) returns contents of matDensity as a double

% --- Executes during object creation, after setting all properties.
function matDensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to matDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SFstruct_Callback(hObject, eventdata, handles)
% hObject    handle to SFstruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SFstruct as text
%        str2double(get(hObject,'String')) returns contents of SFstruct as a double


% --- Executes during object creation, after setting all properties.
function SFstruct_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SFstruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function STmin_Callback(hObject, eventdata, handles)
% hObject    handle to STmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of STmin as text
%        str2double(get(hObject,'String')) returns contents of STmin as a double

% --- Executes during object creation, after setting all properties.
function STmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to STmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function STdel_Callback(hObject, eventdata, handles)
% hObject    handle to STdel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of STdel as text
%        str2double(get(hObject,'String')) returns contents of STdel as a double

% --- Executes during object creation, after setting all properties.
function STdel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to STdel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ParetoFraction_Callback(hObject, eventdata, handles)
% hObject    handle to ParetoFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ParetoFraction as text
%        str2double(get(hObject,'String')) returns contents of ParetoFraction as a double


% --- Executes during object creation, after setting all properties.
function ParetoFraction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ParetoFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function Hidden_FlowDist_Filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radCP1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in FlowFile_Button.
function FlowFile_Button_Callback(hObject, eventdata, handles)
% hObject    handle to FlowFile_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PathFilter1 = [pwd '\Input_Files\Flow_Data\*.dat'];
[Flow_FileName,Flow_PathName] = uigetfile(PathFilter1,'Select the ".dat" file containing flow distribution data.','MultiSelect','off');
 
    if ischar(Flow_FileName)
        OpenFile_Error = 0;
        set(handles.FlowDist_Filename,'String',Flow_FileName);
        set(handles.Hidden_FlowDist_Filename,'String',[pwd '\Input_Files\Flow_Data\' Flow_FileName]);
    else
        OpenFile_Error = 1;
        set(handles.FlowDist_Filename,'String','no file selected');
    end
UpdateEstimates(handles) %calls the function to estimate AEP and CF
    
% --- Executes on button press in Save_GUI.
function Save_GUI_Callback(hObject, eventdata, handles)
% hObject    handle to Save_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%allow the user to specify where to save the settings file
PathFilter2 = [pwd '\Input_Files\GUI_Presets\*.fig'];
[GUI_FileName,GUI_PathName] = uiputfile(PathFilter2,'Save your GUI settings');
    if ~ischar(GUI_FileName)
        return; %if user canceled exit this callback
    end
    
%construct the path name of the save location
saveGUIDataName = fullfile(GUI_PathName,GUI_FileName); 
 
%saves the gui data
hgsave(saveGUIDataName);

% --- Executes on button press in Load_GUI.
function Load_GUI_Callback(hObject, eventdata, handles)
% hObject    handle to Load_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%allow the user to choose which settings to load
PathFilter2 = [pwd '\Input_Files\GUI_Presets\*.fig'];
[LoadGUI_filename, LoadGUI_pathname] = uigetfile(PathFilter2, 'Choose the GUI settings file to load');
    if ~ischar(LoadGUI_filename)
        return; %if user canceled exit this callback
    end 
%construct the path name of the file to be loaded
loadDataName = fullfile(LoadGUI_pathname,LoadGUI_filename);

%this is the gui that will be closed once we load the new settings
theCurrentGUI = gcf;  
%load the settings, which creates a new gui
hgload(loadDataName); 
%closes the old gui
close(theCurrentGUI);

function ElementSpacing_SelectionChangeFcn(hObject, eventdata)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'equal_spacing'
        handles.Elm_spacing = '0';
    case 'cosine_spacing'
        handles.Elm_spacing = '1';
    otherwise       
end
%updates the handles structure
guidata(hObject,handles)



function minD_Callback(hObject, eventdata, handles)
% hObject    handle to minD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minD as text
%        str2double(get(hObject,'String')) returns contents of minD as a double


% --- Executes during object creation, after setting all properties.
function minD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxD_Callback(hObject, eventdata, handles)
% hObject    handle to maxD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxD as text
%        str2double(get(hObject,'String')) returns contents of maxD as a double


% --- Executes during object creation, after setting all properties.
function maxD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minP_Callback(hObject, eventdata, handles)
% hObject    handle to minP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minP as text
%        str2double(get(hObject,'String')) returns contents of minP as a double


% --- Executes during object creation, after setting all properties.
function minP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxP_Callback(hObject, eventdata, handles)
% hObject    handle to maxP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxP as text
%        str2double(get(hObject,'String')) returns contents of maxP as a double


% --- Executes during object creation, after setting all properties.
function maxP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function AEPest_Callback(hObject, eventdata, handles)
% hObject    handle to AEPest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AEPest as text
%        str2double(get(hObject,'String')) returns contents of AEPest as a double


% --- Executes during object creation, after setting all properties.
function AEPest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AEPest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CFest_Callback(hObject, eventdata, handles)
% hObject    handle to CFest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CFest as text
%        str2double(get(hObject,'String')) returns contents of CFest as a double


% --- Executes during object creation, after setting all properties.
function CFest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CFest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in PlotContours.
function PlotContours_Callback(hObject, eventdata, handles)
% hObject    handle to PlotContours (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%read in the needed variables from the GUI
SpdSt = str2num(get(handles.SpdSt,'String'));
SpdEnd = str2num(get(handles.SpdEnd,'String'));
SpdDel = str2num(get(handles.SpdDel,'String'));
U = (SpdSt:SpdDel:SpdEnd)';
rho = str2num(get(handles.Rho,'String'));
minD = str2num(get(handles.minD,'String'));
maxD = str2num(get(handles.maxD,'String'));
minP = str2num(get(handles.minP,'String'));
maxP = str2num(get(handles.maxP,'String'));
ProbDist = str2num(handles.ProbDist);
    if ProbDist == 1 %Rayleigh
        U_mean = str2num(get(handles.U_box,'String'));
        pU = (pi./2).*(U./U_mean.^2).*exp(-(pi./4).*(U./U_mean).^2);
        titleStr = ['Rayleigh U = ' num2str(U_mean,'%3.2f') ' m/s'];
    elseif ProbDist == 2 %Weibull
        Weib_k = str2num(get(handles.k_box,'String'));
        Weib_c = str2num(get(handles.c_box,'String'));
        Weib_Umean = str2num(get(handles.Weib_U_box,'String'));
        pU = (Weib_k./Weib_c).*((U./Weib_c).^(Weib_k-1)).*exp(-(U./Weib_c).^Weib_k);
        titleStr = ['Weibull k = ' num2str(Weib_k,'%3.2f') ', c = ' num2str(Weib_c,'%3.2f') ', U = ' num2str(Weib_Umean,'%3.2f') ' m/s'];
    elseif ProbDist == 3 %user defined
        %User Defined Distribution: read in data from file
        fid = fopen(get(handles.Hidden_FlowDist_Filename,'String'),'rt');
        user_data = textscan(fid,'%f %f','HeaderLines',13);
        fclose(fid);
        U_custom = user_data{1};
        pU_custom = user_data{2};
        %Now interpolate the custom p(U) distribution
        pU = interp1(U_custom,pU_custom,U,'pchip');
        %if the spline caused any negative values change them to zero
        pU(pU<0) = 0;
        titleStr = ['Flow data defined in "' get(handles.FlowDist_Filename,'String') '"'];
    end

%Calculate the AEP and CF contours
    D = linspace(minD,maxD,100); %range of Diameters
    P = linspace(minP,maxP,100); %range of Rated Powers
    estCp = 0.45; %assume a constant efficiency over all flow speeds

    estAEP = zeros(length(D),length(P));
    estCF = zeros(length(D),length(P));
    for n = 1:length(D)
       for m = 1:length(P)
       estPower = 0.5.*rho.*0.25*pi*D(n).^2.*U.^3.*estCp./1000; %kW
       estPower(estPower>P(m)) = P(m); %change power curve to maintain constant rated power

       estAEP(n,m) = 8760*SimpInt(U,estPower.*pU,1,length(U));
       estCF(n,m) = 100*estAEP(n,m)/(8760*P(m)); 
       end
    end

%plot the contours
    figure(4)
    cla reset
    set(gcf,'color','white')
    
    [X Y] = meshgrid(D,P);
    
    hold on
    %plot the CF contours
    valMin = min(min(estCF)); valMin = valMin(1);
    valMax = max(max(estCF)); valMax = valMax(1);
    values = hcosspace(valMin,valMax,20,0);  
    
    [C,h] = contour(X,Y,estCF',values(2:end),'-r');
        handle=clabel(C,h,'fontsize',11,'LabelSpacing',72*4);
        for a=1:length(handle)
            s = get(handle(a),'String'); % get string
            s = str2num(s); % convert in to number
            s = sprintf('%2.1f',s); % format as you need
            set(handle(a),'String',s,'Color','r'); % place it back in the figure
        end

    %plot the AEP contours
    valMin = min(min(estAEP)); valMin = valMin(1);
    valMax = max(max(estAEP)); valMax = valMax(1);
    values = hcosspace(valMin,valMax,20,0);

    [C,h] = contour(X,Y,estAEP',values(2:end),'-k');
        handle=clabel(C,h,'fontsize',11,'LabelSpacing',72*4);
        for a=1:length(handle)
            s = get(handle(a),'String'); % get string
            s = str2num(s); % convert in to number
            s = sprintf('%3.0f',s); % format as you need
            set(handle(a),'String',s,'Color','k'); % place it back in the figure
        end
    hold off
    legend('Capacity Factor (%)','AEP (kWh/yr)',2)
    xlabel('Rotor Diameter (m)');
    ylabel('Rated Power (kW)');
    title(['Estimated contours: ' titleStr],'interpreter','none')
    box on

function UpdateEstimates(handles)
    %read in the needed variables from the GUI
    SpdSt = str2num(get(handles.SpdSt,'String'));
    SpdEnd = str2num(get(handles.SpdEnd,'String'));
    SpdDel = str2num(get(handles.SpdDel,'String'));
    U = (SpdSt:SpdDel:SpdEnd)'; %range of flow speeds
    rho = str2num(get(handles.Rho,'String'));
    Dia = str2num(get(handles.RotorDia,'String'));
    Prated = str2num(get(handles.RatedPwr,'String'));
    ProbDist = str2num(handles.ProbDist);
        U_mean = str2num(get(handles.U_box,'String'));
        Weib_k = str2num(get(handles.k_box,'String'));
        Weib_c = str2num(get(handles.c_box,'String'));
        Weib_Umean = str2num(get(handles.Weib_U_box,'String'));
    opened_file = get(handles.FlowDist_Filename,'String');
    OptMethod = str2num(handles.OptMethod);
    
    if OptMethod == 0 && ProbDist == 4 %opt efficiency is selected
            set(handles.AEPest,'String','n/a');
            set(handles.CFest,'String','n/a');
            return    
    else
        if ProbDist == 1 && U_mean > 0 %Rayleigh
            pU = (pi./2).*(U./U_mean.^2).*exp(-(pi./4).*(U./U_mean).^2);
        elseif ProbDist == 2 && Weib_k > 0 && Weib_c > 0 %Weibull
            pU = (Weib_k./Weib_c).*((U./Weib_c).^(Weib_k-1)).*exp(-(U./Weib_c).^Weib_k);
        elseif ProbDist == 3 && ~strcmp(opened_file,'no file selected')
            %User Defined Distribution: read in data from file
            fid = fopen(get(handles.Hidden_FlowDist_Filename,'String'),'rt');
            user_data = textscan(fid,'%f %f','HeaderLines',13);
            fclose(fid);
            U_custom = user_data{1};
            pU_custom = user_data{2};
            %Now interpolate the custom p(U) distribution
            pU = interp1(U_custom,pU_custom,U,'pchip');
            %if the spline caused any negative values change them to zero
            pU(pU<0) = 0;
        elseif ProbDist == 3 && strcmp(opened_file,'no file selected')
            set(handles.AEPest,'String','n/a');
            set(handles.CFest,'String','n/a');
            return            
        elseif ProbDist == 4 %nothing selected
            set(handles.AEPest,'String','n/a');
            set(handles.CFest,'String','n/a');
            return
        else
            set(handles.AEPest,'String','n/a');
            set(handles.CFest,'String','n/a');
            return            
        end
    end
    
    %Check for any other errors in any of the required variables
    UserError = 0; %set error initially to 0
    if SpdSt<0 || SpdSt>=SpdEnd; UserError = 1; end
    if ProbDist == 3 && ~strcmp(opened_file,'no file selected') && SpdSt<U_custom(1); UserError = 1; end
    if ProbDist == 3 && ~strcmp(opened_file,'no file selected') && SpdEnd>U_custom(end); UserError = 1; end    
    int = (SpdEnd - SpdSt)/SpdDel;
    if SpdDel <= 0;
        UserError = 1;
    elseif SpdDel >= (SpdEnd - SpdSt);
        UserError = 1;
    elseif abs(round(int) - int) > 1e-6; %Checks if integer
        UserError = 1;
    end
    if rho <= 0; UserError = 1; end
    if Dia <= 0; UserError = 1; end
    if Prated <= 0; UserError = 1; end
     
    if UserError == 1
        set(handles.AEPest,'String','n/a');
        set(handles.CFest,'String','n/a');
        return
    else
    %estimate the AEP and CF
       estCp = 0.45; %assume a constant efficiency over all flow speeds
       estPower = 0.5.*rho.*0.25*pi*Dia.^2.*U.^3.*estCp./1000; %kW
       estPower(estPower>Prated) = Prated; %change power curve to maintain constant rated power
       estAEP = 8760*SimpInt(U,estPower.*pU,1,length(U));
       estCF = 100*estAEP/(8760*Prated); 
    %update the text boxes
       AEPstr = num2str(estAEP,'%2.2e');
       newAEPstr = str2num([AEPstr(1) AEPstr(2) AEPstr(3) AEPstr(4)])*10^str2num(AEPstr(end));
       newAEPstr = ThousandSep(newAEPstr);
       
       set(handles.AEPest,'String',newAEPstr);
       set(handles.CFest,'String',num2str(estCF,'%3.1f'));
    end
     
% --- Executes on button press in RecFail.
function RecFail_Callback(hObject, eventdata, handles)
% hObject    handle to RecFail (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


