from basemodel import BaseModel
import siepy
import alphapy


class Alpha(BaseModel):
    def __init__(self, b, x0, y0, e, te, s, alpha):
        super(Alpha, self).__init__(x0, y0, e, te)
        self.b = b
        self.alpha = alpha
        # replaces core radius from s==0 -> 1e-4, fixes /0 situations in potential calculation.
        self.s = s if s != 0.0 else 1e-4

    def modelargs(self, alpha=False):
        if alpha:
            return [self.b, self.x0, self.y0, self.e, self.te, self.s, self.alpha]
        else:
            return [self.b, self.x0, self.y0, self.e, self.te, self.s]

    @BaseModel.standard_frame_rotation
    def phiarray(self, x, y, numexpr=True, *args, **kwargs):
        modelargs = self.modelargs()
        modelargs_with_alpha = self.modelargs(alpha=True)

        if self.alpha == 1.0:
            if self.e == 0.0:
                return siepy.spherical(x, y, modelargs, numexpr=numexpr)
            else:
                return siepy.elliptical(x, y, modelargs, numexpr=numexpr)
        elif self.alpha == -1.0:
            return alphapy.plummer(x, y, modelargs)
        else:
            try:
                from fortran.alpha import alphaf
                return alphaf.general(x, y, modelargs_with_alpha)
            except:
                raise Exception("Alpha!=(0 | -1) not implemented yet")
