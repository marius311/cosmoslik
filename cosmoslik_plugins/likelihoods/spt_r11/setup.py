if __name__=='__main__':

    import os
    
    if not os.path.exists('spt_multif_0809'):
        os.system("wget http://pole.uchicago.edu/public/data/reichardt11/bandpowers_spt_multif_0809.tar.gz")
        os.system("tar zxvf bandpowers_spt_multif_0809.tar.gz")

