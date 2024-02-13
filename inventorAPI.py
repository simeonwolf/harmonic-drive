from win32com.client import constants, Dispatch, GetActiveObject, gencache, CastTo
import numpy as np
from tqdm import tqdm

import matplotlib.pyplot as plt

class InventorAPI():
    '''
    This class adds a curve (xy-Array) to the sketch of an Inventor-Part
    
    In an open Inventor application ist creates a new part, creates a sketch in
    the x-y plane and adds a curve/points to the sketch. The points of the curve
    are transferred to the InventorAPI.pointsTo[Points/Spline/Line](pnts)-function.
    The points must be in the form of a ([x,y],)-vector in an (:,2)-array. The
    points can transfered as points, or connected ba a spline or by a line(--> polygon).
    
    Unit: x,y in [mm]
    
    ! The template parameter must be customized. It has to be the path to the
    template part of your Inventor application.
    
    ! Connecting too many points (~100) to a spline can lead to long calculation
    times in Inventor.
    '''
    def __init__(self):
        try:
            self.oApp = GetActiveObject('Inventor.Application')
        except:
            self.oApp = Dispatch('Inventor.Application')
            self.oApp.Visible = True 
        gencache.EnsureModule('{D98A091D-3A0F-4C3E-B36E-61F62068D488}', 0, 1, 0)
        
        #to be modified:
        template = "C:\\Users\\Public\\Documents\\Autodesk\\Inventor 2024\\Templates\\de-DE\Standard.ipt"
        
        self.invDoc = self.oApp.Documents.Add(constants.kPartDocumentObject, template, True)
        self.tg = self.oApp.TransientGeometry
        self.invPartDoc = CastTo(self.invDoc, 'PartDocument')
        
        xy_plane = self.invPartDoc.ComponentDefinition.WorkPlanes.Item(3)
        self.sketch = self.invPartDoc.ComponentDefinition.Sketches.Add(xy_plane)
        
        print("---Connected to Inventor---")
        
    def pointsToPoints(self, pnts):
        print("---Start transferring the points to the Inventor sketch---")

        for idx in tqdm(range(len(pnts))):
            
            x = pnts[idx, 0]/10
            y = pnts[idx, 1]/10

            self.sketch.SketchPoints.Add( self.tg.CreatePoint2d(x,y) )
            
    def pointsToSpline(self, pnts):
        print("---Start transferring the points to the Inventor sketch---")

        pnts_inventor = self.oApp.TransientObjects.CreateObjectCollection()
        for idx in tqdm(range(len(pnts))):
            
            x = pnts[idx, 0]/10
            y = pnts[idx, 1]/10

            pnts_inventor.Add( self.tg.CreatePoint2d(x,y) )
            
        print("---Creating the spline---")
        self.sketch.SketchSplines.Add(pnts_inventor, 26371)
    
    def pointsToLine(self, pnts, closed = False):
        print("---Start transferring the points to the Inventor sketch---")

        if closed:
            n = len(pnts)+1
        else:
            n = len(pnts)

        for idx in tqdm(range(1,n)):
            
            x = pnts[idx-1, 0]/10
            y = pnts[idx-1, 1]/10
            p_1 = self.tg.CreatePoint2d(x,y)
            
            if closed:
                if idx == len(pnts): idx = 0
            
            x = pnts[idx, 0]/10
            y = pnts[idx, 1]/10
            p_2 = self.tg.CreatePoint2d(x,y)
            
            self.sketch.SketchLines.AddByTwoPoints(p_1, p_2)
                        
def main():
    #test run
    api = InventorAPI()
    
    T = np.linspace(0,2*np.pi,100,endpoint=False)
    
    pnts = np.array([[2*(1-np.cos(t))*np.cos(t),
                      2*(1-np.cos(t))*np.sin(t)] for t in T])
    
    #Pyplot:
    ax = plt.gca()
    ax.axis('equal')
    
    ax.plot(pnts[:,0], pnts[:,1])
    
    #Transfer to Inventor
    api.pointsToSpline(pnts)
    #api.pointsToLine(pnts, closed=True)
    #api.pointsToPoints(pnts)
    
if __name__ == "__main__":
    main()