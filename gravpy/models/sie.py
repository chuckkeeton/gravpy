from basemodel import BaseModel
import siepy


class SIE(BaseModel):
    def __init__(self, b, x0, y0, e, te, s):
        super(SIE, self).__init__(x0, y0, e, te)
        self.b = b
        # replaces core radius from s==0 -> 1e-4, fixes /0 situations in potential calculation.
        self.s = s if s != 0.0 else 1e-4

    def modelargs(self):
        return [self.b, self.x0, self.y0, self.e, self.te, self.s]

    @BaseModel.standard_frame_rotation
    def phiarray(self, x, y, numexpr=True, *args, **kwargs):
        modelargs = self.modelargs()

        if self.e == 0:
            return siepy.spherical(x, y, modelargs, numexpr=numexpr)
        else:
            return siepy.elliptical(x, y, modelargs, numexpr=numexpr)