#ifndef __PHOTONCONTAINERCLONE_H__
#define __PHOTONCONTAINERCLONE_H__

#include <PhotonContainer.h>

class PhotonContainerClone: public PhotonContainer
{
  public:
    PhotonContainerClone(PhotonContainer *photoncont);
    virtual ~PhotonContainerClone() {}
};

#endif /* __PHOTONCONTAINERCLONE_H__ */
