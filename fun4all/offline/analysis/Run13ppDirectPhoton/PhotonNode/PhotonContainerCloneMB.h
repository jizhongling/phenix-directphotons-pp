#ifndef __PHOTONCONTAINERCLONEMB_H__
#define __PHOTONCONTAINERCLONEMB_H__

#include <PhotonContainerMB.h>

class PhotonContainerCloneMB: public PhotonContainerMB
{
  public:
    PhotonContainerCloneMB(PhotonContainerMB *photoncont);
    virtual ~PhotonContainerCloneMB() {}
};

#endif /* __PHOTONCONTAINERCLONEMB_H__ */
