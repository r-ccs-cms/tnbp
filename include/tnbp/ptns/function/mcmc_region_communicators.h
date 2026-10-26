/// This file is a part of r-ccs-cms/tnbp
#ifndef TNBP_PTNS_FUNCTION_MCMC_REGION_COMMUNICATORS_H
#define TNBP_PTNS_FUNCTION_MCMC_REGION_COMMUNICATORS_H
#include <mpi.h>
#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
namespace tnbp {
/** The caller supplies the same sorted region sets on every parent rank.
 * Only members call create_group. No parent-wide collective is performed.
 * Handles are owned by the sampling call and are reused across all sweeps.
 */
inline void mcmc_acquire_region_communicators(MPI_Comm parent,int creation_tag,
    const std::vector<std::vector<int>>& sets,std::map<std::vector<int>,MPI_Comm>& handles) {
    int rank,size;MPI_Comm_rank(parent,&rank);MPI_Comm_size(parent,&size);
    int present=0;int* upper=nullptr;MPI_Comm_get_attr(parent,MPI_TAG_UB,&upper,&present);
    if(creation_tag<0 || !present || creation_tag>*upper)
        throw std::invalid_argument("invalid region communicator creation tag");
    if(!std::is_sorted(sets.begin(),sets.end()) || std::adjacent_find(sets.begin(),sets.end())!=sets.end())
        throw std::invalid_argument("region owner sets must be sorted and unique");
    // Validate before the first communication, including nonparticipant inputs.
    for(const auto& members:sets) {
        if(members.empty() || !std::is_sorted(members.begin(),members.end()) ||
           std::adjacent_find(members.begin(),members.end())!=members.end() || members.front()<0 || members.back()>=size)
            throw std::invalid_argument("invalid region owner set");
    }
    for(const auto& members:sets) {
        if(handles.count(members))continue;
        MPI_Comm local=MPI_COMM_NULL;
        if(std::binary_search(members.begin(),members.end(),rank)) {
            MPI_Group group,subset;MPI_Comm_group(parent,&group);
            MPI_Group_incl(group,static_cast<int>(members.size()),members.data(),&subset);
            MPI_Comm_create_group(parent,subset,creation_tag,&local);
            MPI_Group_free(&subset);MPI_Group_free(&group);
        }
        handles.emplace(members,local);
    }
}
inline void mcmc_release_region_communicators(std::map<std::vector<int>,MPI_Comm>& handles) {
    for(auto& entry:handles)if(entry.second!=MPI_COMM_NULL)MPI_Comm_free(&entry.second);
    handles.clear();
}
template<class Site>
std::map<Site,Site> mcmc_region_local_owners(const std::map<Site,Site>& original,
    const std::vector<Site>& sites,const std::vector<int>& members) {
    std::map<Site,Site> result;
    for(auto site:sites) {
        const auto rank=original.at(site);
        const auto it=std::lower_bound(members.begin(),members.end(),rank);
        if(it==members.end() || *it!=rank)throw std::invalid_argument("missing region owner");
        result.emplace(site,static_cast<Site>(it-members.begin()));
    }
    return result;
}
} // namespace tnbp
#endif
