import unittest

import Region1basic as bas1
import Region1backward as bk1

class Region1Test (unittest.TestCase):

      def setUp(self):
         # IF97-Rev,Table 5 Page 9 : p,T(K),h(kJ/kg),v(m^3/kg),u(kJ/kg),s(kJ/(kg.k)),Cp(kJ/(kg.k)),w(m/s)
         self.tab5=[ [  3,  300,  0.115331273e3 ,0.100215168e-2,  0.112324818e3, 0.392294792,  0.417301218e1 ,  0.150773921e4],
                   [80,  300,  0.184142828e3 , 0.971180894e-3 ,  0.106448356e3,  0.368563852,  0.401008987e1,  0.163469054e4],
                   [ 3,  500,  0.975542239e3 , 0.120241800e-2 ,  0.971934985e3,  0.258041912e1,  0.465580682e1,  0.124071337e4]]
 
         
        # IF97-Rev, Table 7.Page 11 p( MPa),h(kJ/kg), T(K)
         self.tab7=[ [  3,  500,  0.391798509e3 ],
                   [80,  500,    0.378108626e3 ],
                   [ 80,   1500,  0.611041229e3]]

         # IF97-Rev, Table 9.Page 11 p( MPa),s(kJ/(kg.k), T(K)
         self.tab9=[ [  3,  0.5, 0.307842258e3 ],
                   [80,  0.5,  0.309979785e3 ],
                   [80,  3,  0.565899909e3]]

         # Supp-PHS12-2014, Table 3, Page 6 h(kJ/kg),s(kJ/(kg.k), p( MPa)
         self.PHS12_tab3=[ [0.001,  0, 9.800980612e-4],
                           [90,     0, 9.192954727e1 ],
                          [1500,   3.4,  5.868294423e1]]
         

   
      def testSpecificEnthalpy_PT(self):
           places = 6
           for item in self.tab5:
               self.assertAlmostEqual(bas1.enthalpyreg1(item[0], item[1]),item[2],places)

      def testSpecificVolume_PT(self):
           places = 6
           for item in self.tab5:
               self.assertAlmostEqual(bas1.volreg1(item[0], item[1]),item[3], places)
               
      def testSpecificinternalEnergy_PT(self):
           places = 6
           for item in self.tab5:
               self.assertAlmostEqual(bas1.energyreg1(item[0], item[1]),item[4],places) 
               
      def testSpecificEntropy_PT(self):
           places = 6
           for item in self.tab5:
               self.assertAlmostEqual(bas1.entropyreg1(item[0], item[1]),item[5],places) 
               
      def testSpecificisobaricheatcapacity_PT(self):
           places = 6
           for item in self.tab5:
               self.assertAlmostEqual(bas1.cpreg1(item[0], item[1]),item[6],places)  
               
      def testSpeedofsound_PT(self):
           places = 5
           for item in self.tab5:
               self.assertAlmostEqual(bas1.spsoundreg1(item[0], item[1]),item[7],places)                                      
    
      def testTemperature_ph(self):
           places = 6
           for item in self.tab7:
               self.assertAlmostEqual(bk1.phtoTreg1(item[0], item[1]),item[2],places)
             
      def testTemperature_ps(self):
           places = 6
           for item in self.tab9:
               self.assertAlmostEqual(bk1.pstoTreg1(item[0], item[1]),item[2],places)
             
      def testPressure_hs(self):
           places = 8
           for item in self.PHS12_tab3:
               self.assertAlmostEqual(bk1.hstopreg1(item[0], item[1]),item[2],places)
             

def suite_use_make_suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Region1Test))
    return suite

if __name__ == '__main__':
    unittest.main()  
