import win32com.client
from win32com.client import constants
import numpy as np

class InventorAPI():
    def __init__(self):
        self.oApp = win32com.client.Dispatch('Inventor.Application')
        
        template = "C:\\Users\\Public\\Documents\\Autodesk\\Inventor 2024\\Templates\\de-DE\Standard.ipt"
        
        #template = "C:\Autodesk\Inventor_2022" + chr(92) + "templates\igus_Bauteil.ipt"
        self.invDoc = self.oApp.Documents.Add(constants.kPartDocumentObject, template, True)
        self.tg = self.oApp.TransientGeometry
        self.invPartDoc = win32com.client.CastTo(self.invDoc, 'PartDocument')
        print("---Verbunden mit Inventor---")
        
    def pointsToSpline(self, pnts):
        print("---Start der Übetragung der Punkte in Inventor-Skizze---")
        
        xy_plane = self.invPartDoc.ComponentDefinition.WorkPlanes.Item(3)
        sketch = self.invPartDoc.ComponentDefinition.Sketches.Add(xy_plane)
        
        pnts_inventor = self.oApp.TransientObjects.CreateObjectCollection()

        for idx in range(len(pnts)):
            
            x = pnts[idx, 0]/10
            y = pnts[idx, 1]/10
            
            print("---", idx/len(pnts)*100, "% ---")

            #sketch.SketchPoints.Add(self.tg.CreatePoint2d(x,y))

            pnts_inventor.Add( self.tg.CreatePoint2d(x,y) )
         
        print("---Erzeugung des Splines---")
            
        sketch.SketchSplines.Add(pnts_inventor, 26371)

def main():
    inventor_api = InventorAPI()
    
    punkte = np.array([[1,2], [3,4]])
    
    inventor_api.pointsToSpline(punkte)
    
if __name__ == "__main__":
    main()

