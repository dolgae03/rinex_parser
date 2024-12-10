function [Obs_columns, nObs_types] = obs_type_find(Obs_types, sysId)

% SYNTAX:
%   [Obs_columns, nObs_types] = obs_type_find(Obs_types, sysId);
%
% INPUT:
%   Obs_types = cell of strings containing observation types
%               RINEX v2.xx --> e.g. L1C1P1...
%               RINEX v3.xx --> e.g. C1CL1CD1C...
%   sysId = cell-array containing one-letter identifiers for constellations
%
% OUTPUT:
%   Obs_columns = structure containing the column number of each observation type
%                 in the following fields:
%                   .L1 = L1 column
%                   .L2 = L2 column
%                   .C1 = C1 column
%                   .P1 = P1 column
%                   .P2 = P2 column
%                   .S1 = S1 column
%                   .S2 = S2 column
%                   .D1 = D1 column
%                   .D2 = D2 column
%                 In the case of RINEX v3.xx, an additional field is added
%                 for specifying the constellation, e.g.:
%                   .G.L1 (GPS)
%                   .R.L1 (GLONASS)
%   nObs_types = number of available observation types
%
% DESCRIPTION:
%   Detection of the column index for phase observations (L1, L2), for
%   code observations (C1, P1, P2), SNR ratios (S1, S2) and Doppler
%   measurements (D1, D2).

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Stefano Caldera
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

if (isempty(sysId)) %RINEX v2.xx

    nObs_types = size(Obs_types{1},2)/2;

    %search L1 column
    s1 = strfind(Obs_types{1}, 'L1');
    s2 = strfind(Obs_types{1}, 'LA');
    s = [s1 s2];
    col_L1 = (s+1)/2;

    %search L2 column
    s1 = strfind(Obs_types{1}, 'L2');
    s2 = strfind(Obs_types{1}, 'LC');
    s = [s1 s2];
    col_L2 = (s+1)/2;

    %search C1 column
    s1 = strfind(Obs_types{1}, 'C1');
    s2 = strfind(Obs_types{1}, 'CA');
    s = [s1 s2];
    col_C1 = (s+1)/2;

    %search P1 column
    s1 = strfind(Obs_types{1}, 'P1');
    s2 = strfind(Obs_types{1}, 'CA'); %QZSS does not use P1
    s = [s1 s2];
    col_P1 = (s+1)/2;

    %if RINEX v2.12 and GPS/GLONASS P1 observations are not available
    if (length(col_P1) ~= 2 && ~isempty(s2))
        %keep QZSS CA observations as C1
        col_P1 = [];
    end

    %search P2 column
    s1 = strfind(Obs_types{1}, 'P2');
    s2 = strfind(Obs_types{1}, 'CC');
    s = [s1 s2];
    col_P2 = (s+1)/2;

    %search S1 column
    s1 = strfind(Obs_types{1}, 'S1');
    s2 = strfind(Obs_types{1}, 'SA');
    s = [s1 s2];
    col_S1 = (s+1)/2;

    %search S2 column
    s1 = strfind(Obs_types{1}, 'S2');
    s2 = strfind(Obs_types{1}, 'SC');
    s = [s1 s2];
    col_S2 = (s+1)/2;

    %search D1 column
    s1 = strfind(Obs_types{1}, 'D1');
    s2 = strfind(Obs_types{1}, 'DA');
    s = [s1 s2];
    col_D1 = (s+1)/2;

    %search D2 column
    s1 = strfind(Obs_types{1}, 'D2');
    s2 = strfind(Obs_types{1}, 'DC');
    s = [s1 s2];
    col_D2 = (s+1)/2;

    Obs_columns.L1 = col_L1;
    Obs_columns.L2 = col_L2;
    Obs_columns.C1 = col_C1;
    Obs_columns.P1 = col_P1;
    Obs_columns.P2 = col_P2;
    Obs_columns.S1 = col_S1;
    Obs_columns.S2 = col_S2;
    Obs_columns.D1 = col_D1;
    Obs_columns.D2 = col_D2;

else %RINEX v3.xx
    for c = 1 : length(sysId)

        nObs_types.(sysId{c}) = size(Obs_types.(sysId{c}),2)/3;

        switch sysId{c}
            case 'G' %GPS
                idC1 = {'C1C'};               %L1
                idP1 = {'C1W';'C1P'};         %L1
                idL1 = {'L1W';'L1X';'L1C'};   %L1
                idS1 = {'S1W';'S1P';'S1C'};   %L1
                idD1 = {'D1W';'D1P';'D1C'};   %L1
                %--------------------------------
                %idC2 = {'C2W';'C2P';'C2C};    %L2,L2C (precedence to L2)
                idP2 = {'C2W';'C2P';'C2C';'C2S';'C2L';'C2X'};    %L2,L2C
                idL2 = {'L2W';'L2P';'L2C';'L2S';'L2L';'L2X'};    %L2,L2C
                idS2 = {'S2W';'S2P';'S2C';'S2S';'S2L';'S2X'};    %L2,L2C
                idD2 = {'D2W';'D2P';'D2C';'D2S';'D2L';'D2X'};    %L2,L2C
                %--------------------------------
                idP3 = {'C5Q'};               %L5
                idL3 = {'L5Q'};               %L5
                idS3 = {'S5Q'};               %L5
                idD3 = {'D5Q'};               %L5
                %--------------------------------
                idP4 = {};               
                idL4 = {};               
                idS4 = {};               
                idD4 = {};                               
                %--------------------------------
                idP5 = {};               
                idL5 = {};               
                idS5 = {};               
                idD5 = {};                               
                %--------------------------------
            case 'R' %GLONASS
                idC1 = {'C1C'};               %R1
                idP1 = {'C1P'};               %R1
                idL1 = {'L1C'};               %R1
                idS1 = {'S1C'};               %R1
                idD1 = {'D1C'};               %R1
                %--------------------------------
                idP2 = {'C2P'};               %R2
                idL2 = {'L2P'};               %R2
                idS2 = {'S2P'};               %R2
                idD2 = {'D2P'};               %R2
                %--------------------------------
                idP3 = {};               %L5
                idL3 = {};               %L5
                idS3 = {};               %L5
                idD3 = {};               %L5
                %--------------------------------
                idP4 = {};               
                idL4 = {};               
                idS4 = {};               
                idD4 = {};
                %--------------------------------
                idP5 = {};               
                idL5 = {};               
                idS5 = {};               
                idD5 = {};                               
                %--------------------------------
            case 'E' %Galileo
                idC1 = {'C1X';'C1C'};         %E1
                idP1 = {'...'};               %E1
                idL1 = {'L1X';'L1C'};         %E1
                idS1 = {'S1X';'S1C'};         %E1
                idD1 = {'D1X';'D1C'};         %E1
                %--------------------------------
                %idP2 = {'C5X';'C5Q'};        %E5a
                %idL2 = {'L5X';'L5Q'};        %E5a
                %idS2 = {'S5X';'S5Q'};        %E5a
                %idD2 = {'D5X';'D5Q'};        %E5a
                %--------------------------------
                idP2 = {'C7X';'C7Q'};        %E5b
                idL2 = {'L7X';'L7Q'};        %E5b
                idS2 = {'S7X';'S7Q'};        %E5b
                idD2 = {'D7X';'D7Q'};        %E5b
                %--------------------------------
                %idP2 = {'C8X';'C8Q'};         %E5
                %idL2 = {'L8X';'L8Q'};         %E5
                %idS2 = {'S8X';'S8Q'};         %E5
                %idD2 = {'D8X';'D8Q'};         %E5
                %--------------------------------
                idP4 = {};               
                idL4 = {};               
                idS4 = {};               
                idD4 = {};          
                %--------------------------------
                idP5 = {};               
                idL5 = {};               
                idS5 = {};               
                idD5 = {};                               
                %--------------------------------
            case 'C' %Compass/Beidou
                idC1 = {'C1I';'C1Q';};   %B1
                idP1 = {'...'};               %B1
                idL1 = {'L1I';'L1Q';};   %B1
                idS1 = {'S1I';'S1Q';};   %B1
                idD1 = {'D1I';'D1Q';};   %B1

                % idC1 = {'C1I';'C1Q';'C2I'};   %B1
                % idP1 = {'...'};               %B1
                % idL1 = {'L1I';'L1Q';'L2I'};   %B1
                % idS1 = {'S1I';'S1Q';'S2I'};   %B1
                % idD1 = {'D1I';'D1Q';'D2I'};   %B1
                %--------------------------------
                idL2 = {'L7I';'L7Q'};         %B2b
                idP2 = {'C7I';'C7Q'};         %B2b
                idS2 = {'S7I';'S7Q'};         %B2b
                idD2 = {'D7I';'D7Q'};         %B2b
                %--------------------------------
                idL3 = {'L5D';'L5P'};         %B2a
                idP3 = {'C5D';'C5P'};         %B2a
                idS3 = {'S5D';'S5P'};         %B2a
                idD3 = {'D5D';'D5P'};         %B2a
                %--------------------------------
                idL4 = {'L7D';'L7P'};         %B2b           
                idP4 = {'C7D';'C7P'};         %B2b      
                idS4 = {'S7D';'S7P'};         %B2b      
                idD4 = {'D7D';'D7P'};         %B2b  
                %--------------------------------
                idP5 = {'C1P'};               
                idL5 = {'L1P'};               
                idS5 = {'S1P'};               
                idD5 = {'D1P'};                               
                %--------------------------------
                
            case 'J' %QZSS
                idC1 = {'C1X';'C1C'};         %J1
                idP1 = {'...'};               %J1
                idL1 = {'L1X';'L1C'};         %J1
                idS1 = {'S1X';'S1C'};         %J1
                idD1 = {'D1X';'D1C'};         %J1
                %--------------------------------
                idP2 = {'C2X';'C2C';'C2S'};         %J2
                idL2 = {'L2X';'L2C';'L2S'};         %J2
                idS2 = {'S2X';'S2C';'S2S'};         %J2
                idD2 = {'D2X';'D2C';'D2S'};         %J2
                %--------------------------------
                idP3 = {'C5Q'};         %J5
                idL3 = {'L5Q'};         %J5
                idS3 = {'S5Q'};         %J5
                idD3 = {'D5Q'};         %J5
                %--------------------------------
                idP4 = {};               
                idL4 = {};               
                idS4 = {};               
                idD4 = {};                               
                %--------------------------------
                idP5 = {};               
                idL5 = {};               
                idS5 = {};               
                idD5 = {};                               
                %--------------------------------
        end

        %search L1 column
        s=[];
        for i = 1 : length(idL1)
            s = strfind(Obs_types.(sysId{c}), idL1{i}); if (~isempty(s)), break, end;
        end
        col_L1 = (s+2)/3;

        %search L2 column
        s=[];
        for i = 1 : length(idL2)
            s = strfind(Obs_types.(sysId{c}), idL2{i}); if (~isempty(s)), break, end;
        end
        col_L2 = (s+2)/3;
        
        %search L3 column
        s=[];
        for i = 1 : length(idL3)
            s = strfind(Obs_types.(sysId{c}), idL3{i}); if (~isempty(s)), break, end;
        end
        col_L3 = (s+2)/3; 
        
        %search L4 column
        s=[];
        for i = 1 : length(idL4)
            s = strfind(Obs_types.(sysId{c}), idL4{i}); if (~isempty(s)), break, end;
        end
        col_L4 = (s+2)/3; 

        %search L5 column
        s=[];
        for i = 1 : length(idL5)
            s = strfind(Obs_types.(sysId{c}), idL5{i}); if (~isempty(s)), break, end;
        end
        col_L5 = (s+2)/3; 
        
        %search C1 column
        s=[];
        for i = 1 : length(idC1)
            s = strfind(Obs_types.(sysId{c}), idC1{i}); if (~isempty(s)), break, end;
        end
        col_C1 = (s+2)/3;

        %search P1 column
        s=[];
        for i = 1 : length(idP1)
            s = strfind(Obs_types.(sysId{c}), idP1{i}); if (~isempty(s)), break, end;
        end
        col_P1 = (s+2)/3;

        %search P2 column
        s=[];
        for i = 1 : length(idP2)
            s = strfind(Obs_types.(sysId{c}), idP2{i}); if (~isempty(s)), break, end;
        end
        col_P2 = (s+2)/3;
        
        %search P3 column
        s=[];
        for i = 1 : length(idP3)
            s = strfind(Obs_types.(sysId{c}), idP3{i}); if (~isempty(s)), break, end;
        end
        col_P3 = (s+2)/3;

        %search P4 column
        s=[];
        for i = 1 : length(idP4)
            s = strfind(Obs_types.(sysId{c}), idP4{i}); if (~isempty(s)), break, end;
        end
        col_P4 = (s+2)/3;

        %search P5 column
        s=[];
        for i = 1 : length(idP5)
            s = strfind(Obs_types.(sysId{c}), idP5{i}); if (~isempty(s)), break, end;
        end
        col_P5 = (s+2)/3;
        
        %search S1 column
        s=[];
        for i = 1 : length(idS1)
            s = strfind(Obs_types.(sysId{c}), idS1{i}); if (~isempty(s)), break, end;
        end
        col_S1 = (s+2)/3;

        %search S2 column
        s=[];
        for i = 1 : length(idS2)
            s = strfind(Obs_types.(sysId{c}), idS2{i}); if (~isempty(s)), break, end;
        end
        col_S2 = (s+2)/3;
        
        %search S3 column
        s=[];
        for i = 1 : length(idS3)
            s = strfind(Obs_types.(sysId{c}), idS3{i}); if (~isempty(s)), break, end;
        end
        col_S3 = (s+2)/3;

        %search S4 column
        s=[];
        for i = 1 : length(idS4)
            s = strfind(Obs_types.(sysId{c}), idS4{i}); if (~isempty(s)), break, end;
        end
        col_S4 = (s+2)/3;       

        %search S5 column
        s=[];
        for i = 1 : length(idS5)
            s = strfind(Obs_types.(sysId{c}), idS5{i}); if (~isempty(s)), break, end;
        end
        col_S5 = (s+2)/3;    
        
        %search D1 column
        s=[];
        for i = 1 : length(idD1)
            s = strfind(Obs_types.(sysId{c}), idD1{i}); if (~isempty(s)), break, end;
        end
        col_D1 = (s+2)/3;

        %search D2 column
        s=[];
        for i = 1 : length(idD2)
            s = strfind(Obs_types.(sysId{c}), idD2{i}); if (~isempty(s)), break, end;
        end
        col_D2 = (s+2)/3;
        
        %search D3 column
        s=[];
        for i = 1 : length(idD3)
            s = strfind(Obs_types.(sysId{c}), idD3{i}); if (~isempty(s)), break, end;
        end
        col_D3 = (s+2)/3;      

        %search D4 column
        s=[];
        for i = 1 : length(idD4)
            s = strfind(Obs_types.(sysId{c}), idD4{i}); if (~isempty(s)), break, end;
        end
        col_D4 = (s+2)/3;        

        %search D4 column
        s=[];
        for i = 1 : length(idD5)
            s = strfind(Obs_types.(sysId{c}), idD5{i}); if (~isempty(s)), break, end;
        end
        col_D5 = (s+2)/3;        
        
        Obs_columns.(sysId{c}).L1 = col_L1;
        Obs_columns.(sysId{c}).L2 = col_L2;
        Obs_columns.(sysId{c}).L3 = col_L3;
        Obs_columns.(sysId{c}).L4 = col_L4;
        Obs_columns.(sysId{c}).L5 = col_L5;         
        Obs_columns.(sysId{c}).C1 = col_C1;
        Obs_columns.(sysId{c}).P1 = col_P1;
        Obs_columns.(sysId{c}).P2 = col_P2;
        Obs_columns.(sysId{c}).P3 = col_P3;
        Obs_columns.(sysId{c}).P4 = col_P4;  
        Obs_columns.(sysId{c}).P5 = col_P5;        
        Obs_columns.(sysId{c}).S1 = col_S1;
        Obs_columns.(sysId{c}).S2 = col_S2;
        Obs_columns.(sysId{c}).S3 = col_S3;
        Obs_columns.(sysId{c}).S4 = col_S4;    
        Obs_columns.(sysId{c}).S5 = col_S5;       
        Obs_columns.(sysId{c}).D1 = col_D1;
        Obs_columns.(sysId{c}).D2 = col_D2;
        Obs_columns.(sysId{c}).D3 = col_D3;
        Obs_columns.(sysId{c}).D4 = col_D4;   
        Obs_columns.(sysId{c}).D5 = col_D5;         
    end
end
