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
round(table(bid_aucnum)/668,3)
