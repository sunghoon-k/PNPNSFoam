''' 사용시 주의 사항
경계 조건에 맞는 기초 폴더가 있어야 한다. ex) 'slip_Case', 'Smoluchowski_Case'
생성할 폴더가 있어야 한다. 

높이: 10um에 대해 격자수 158개 ㄱㄱ, grading 1.03 -> thick_max/thick_min = 1.03^157

비율 계산하는 방법: 메쉬 개수 / 직경 제곱

둘레에 부여하는 격자수: dn
정사각 격자일 때 dl = 10e-6/hn/0.9
정사각 격자 면적과 정삼각 격자의 면적을 같도록 => dl^2 = sqrt(3)/4 * dr^2
=> dr = sqrt(4/sqrt(3)) * dl
직경 D에 대해 D*pi = dn * dr
dn = D/dr
'''

# Diameter는 고정하고, Height를 통해 ratio를 조정하자. 왜? 그래야 격자수 최소화 가능
# 가장 큰 경우는 직경 20에 높이 10인 경우 
# 그렇다면 무차원 직경과 높이는? 무차원 높이를 1로 고정하고, 직경은 ratio에 비례하도록 ㄱㄱ
import os
import subprocess
from math import *
import numpy as np

Ratios  = [2, 4, 8, 16] # [1, 2, 4, 8, 16] # 종횡비
h_base = 10 # 단위 [um]
l0 = 20
a1 = 3e-3 # 3e-3 [um] = 3 [nm]

Diameter = 20

r = 1.03 
dl = 1/3 # 1/3 um
dr = np.sqrt(4/np.sqrt(3)) * dl

boundaries = ['slip', 'noSlip'] # ['slip', 'Smoluchowski']

# Height = 5

os.system(f"mkdir Diameter{Diameter}")
# 높이에 부여하는 격자 수는 고정하고 둘레에 부여하는 격자수로 coarse/fine 조정한다. 
Dn = int(Diameter*np.pi/dr/4)
Dn = int(Dn/2) + 1 # coarse: 기존의 1/2 만큼

for Ratio in Ratios: # (Aspect Ratio) 종횡비 
  Height = Diameter / Ratio
  
  os.system(f"mkdir Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}")
  
  Hn = int(log(1 + min(Height, h_base) * (r - 1) / (a1) )/log(r))
  
  str_head = f"""
Ratio = {Ratio};
h_base = {min(Height, 10)};
h = {Height};
Dn1 = {Dn};
hn1 = {Hn}; // number of layer
//rMaxDivideMin = {pow(1.03, Hn-1)};
r = {r}; //rMaxDivideMin ^ (1/(hn1-1)); // layer마다 증가하는 비율
    """
  str_body = str_head + """
/* **************************************************** */

//+
SetFactory("OpenCASCADE");
Point(1) = {-1/2*Ratio, 0         , 0, 1.0};
Point(2) = {0         , 1/2*Ratio , 0, 1.0};
Point(3) = {1/2*Ratio , 0         , 0, 1.0};
Point(4) = {0         , -1/2*Ratio, 0, 1.0};
Point(5) = {-1/3*Ratio, 0         , 0, 1.0};
Point(6) = {0         , 1/3*Ratio , 0, 1.0};
Point(7) = {1/3*Ratio , 0         , 0, 1.0};
Point(8) = {0         , -1/3*Ratio, 0, 1.0};
Point(9) = {0         , 0         , 0, 1.0};

Point(10) = {(Sqrt(2)/2)*1/8*Ratio , -(Sqrt(2)/2)*1/8*Ratio, 0, 1.0};
Point(11) = {-(Sqrt(2)/2)*1/8*Ratio, -(Sqrt(2)/2)*1/8*Ratio, 0, 1.0};
Point(12) = {-(Sqrt(2)/2)*1/8*Ratio, (Sqrt(2)/2)*1/8*Ratio , 0, 1.0};
Point(13) = {(Sqrt(2)/2)*1/8*Ratio , (Sqrt(2)/2)*1/8*Ratio , 0, 1.0};


Circle(1) = {6, 11, 7};
Circle(2) = {7, 12, 8};
Circle(3) = {8, 13, 5};
Circle(4) = {5, 10, 6};
Circle(5) = {2, 9, 3};
Circle(6) = {3, 9, 4};
Circle(7) = {4, 9, 1};
Circle(8) = {1, 9, 2};
Line(9) = {2, 6};
Line(10) = {3, 7};
Line(11) = {4, 8};
Line(12) = {1, 5};
Recursive Delete {
  Point{9, 10, 11, 12, 13}; 
}//+

Line Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};

Line Loop(2) = {1, -10, -5, 9};
Plane Surface(2) = {2};
Line Loop(3) = {2, -11, -6, 10};
Plane Surface(3) = {3};
Line Loop(4) = {3, -12, -7, 11};
Plane Surface(4) = {4};
Line Loop(5) = {4, -9, -8, 12};
Plane Surface(5) = {5};

//+
Transfinite Surface {1} = {5, 6, 7, 8};
//+
Transfinite Surface {5} = {1, 2, 6, 5};
//+
Transfinite Surface {2} = {2, 3, 7, 6};
//+
Transfinite Surface {3} = {3, 4, 8, 7};
//+
Transfinite Surface {4} = {4, 1, 5, 8};
//+
Transfinite Line {4, 1, 2, 3, 12, 9, 10, 11, 8, 5, 6, 7} = Dn1 Using Progression 1;
//+
Recombine Surface {1, 5, 2, 3, 4};

// 2D -> 3D
For i In {0:(hn1-2)} // n-1번 실행
//Printf("%f", 2^i);
	layerList1[i] = 1; // 1장을 깐다. 
	layerList2[i] = (h_base/h) * (r^(i+1) - 1)/(r^hn1 - 1); // Layer의 두께	
EndFor

layerList1[hn1-1] = 1;
layerList2[hn1-1] = h_base/h; // 마지막 layer는 무조건 1

/* 10um 이상 높이 처리하는 부분 */	
If(h>h_base)
a = (h_base/h) * (r-1)/(r^hn1-1); // 초항
a_n = a*r^(hn1-1);
layerList1[hn1] = (1-(h_base/h))/a_n; // 몇 개의 layer로 나누나?
layerList2[hn1] = 1; 
EndIf

	
Extrude {0, 0, 1} {Surface{1, 2, 3, 4, 5}; Layers{layerList1[], layerList2[]}; Recombine;}

// top
Physical Surface("top")    = {10, 14, 17, 20, 22};
// Bottom
Physical Surface("bottom") = {1, 2, 3, 4, 5};
// Wall
Physical Surface("wall")  = {12, 16, 19, 21};

// Physical Volume
Physical Volume("Model")   = {1, 2, 3, 4, 5};    
// GMSH  ----------------------------------- H Version -----------------------------------  GMSH 

// Coherence;

// save  -------------------------------------------------------------------------------------
//Mesh.MshFileVersion = 2.2;
Mesh 3;
    """
  meshTitle = f"Diameter{Diameter}_Ratio{Ratio}_Dn{Dn-1}"
  str_tail = str_body + f"""Save "{meshTitle}.msh";"""

  with open(f'{meshTitle}.geo', 'w') as f:
    f.write(str_tail)
          # f.write(str_body)

  os.system(f"gmsh -nt 20 {meshTitle}.geo -") # {meshTitle}.geo -") # msh파일 제작
  for boundary in boundaries:
    os.system(f"cp -r {boundary}_Case/ Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/{meshTitle}_{boundary}/") # 케이스 폴더 만들기 # ./solverCollection/{title}/
    os.system(f"cp {meshTitle}.msh Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/{meshTitle}_{boundary}/") # msh 파일 옮기기
    
    # alias 사용하기 참조: https: // brownbears.tistory.com / 204
    command1 = f'cd Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/{meshTitle}_{boundary}; fe40; gmshToFoam {meshTitle}.msh' # ; PNPNSFoam_org_Ueof'
    sp = subprocess.Popen(["/bin/bash", "-i", "-c", command1])
    sp.communicate()
  os.system(f"rm {meshTitle}.geo") # geo 파일 삭제
  os.system(f"rm {meshTitle}.msh") # msh 파일 삭제

