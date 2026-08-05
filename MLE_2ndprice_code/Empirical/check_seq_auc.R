raw_XBOX_7=read.csv("Empirical/Xbox_7day_auctions.csv")
bidders_List=unique(raw_XBOX_7$bidder)
bid_aucnum=rep(0,length(bidders_List))
names(bid_aucnum)=bidders_List
bid_list=unique(raw_XBOX_7$auction)
for (bid_id in bid_list) {
  temp=raw_XBOX_7[raw_XBOX_7$auction==bid_id,]
  temp_bidders=unique(temp$bidder)
  for (bidder in temp_bidders) {
    bid_aucnum[bidder]=bid_aucnum[bidder]+1
  }
}
max(bid_aucnum)
mean(bid_aucnum>1)
table(bid_aucnum)
round(table(bid_aucnum)/668,3) ## less than 14%


## check proxy and multiple
bidder_multi=rep(0,length(bid_list))
bidder_multi2=rep(0,length(bid_list))
bid_aucnum_proxy=rep(0,length(bid_list))
for (bid_id in 1:length(bid_list)) {
  temp=raw_XBOX_7[raw_XBOX_7$auction==bid_list[bid_id],]
  temp_bidders=unique(temp$bidder)
  
  count_bidder=rep(0,length(temp_bidders))
  count_bidder_proxy=rep(0,length(temp_bidders))
  
  for (bidder in 1:length(temp_bidders)) {
    bidder_pos=which(temp$bidder==temp_bidders[bidder])
    ## count multiple bids in consecutive 
    if (any(diff(bidder_pos)==1)) {
      count_bidder_proxy[bidder]=1
      ## removing consecutive bids
      bidder_pos=bidder_pos[-which(diff(bidder_pos)==1)]
    }
    ## count multiple bids
    count_bidder[bidder]=length(bidder_pos)
  }
  bidder_multi[bid_id]=sum(count_bidder>1)/length(temp_bidders)
  bidder_multi2[bid_id]=sum(count_bidder>2)/length(temp_bidders)
  bid_aucnum_proxy[bid_id]=sum(count_bidder_proxy)/length(temp_bidders)
}
mean(bidder_multi) ## around 20% after removing consecutive bids
mean(bid_aucnum_proxy) ## less than 5%

mean(bidder_multi2) ## about 6% have more than 2 bids