import dilap.core.context as dgc
import dilap.generate.landscape as dls
import dilap.generate.house as dlh

class lot(dgc.context):

    def generate(self,worn = 0):
        house = dlh.house(35,25)
        #house.generate(worn)
        #uself._consume(house)
        lscape = dls.landscape()
        lscape.generate(worn,house)
        self._consume(lscape)

   
