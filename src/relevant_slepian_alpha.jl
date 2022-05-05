
# Software note: this code was originally developed by F. Simons and D. Wang in
# the matlab programming language, and was rewritten in julia for performance.
# The original code is available at:
# https://github.com/csdms-contrib/slepian_{alpha,foxtrot}
# and is distributed under the GNU GPL v 2.0 license (original license below).
# I take responsibility for any errors in this code as it is not a verbatim copy
# of the original. In particular, some inbuilt matlab functions (e.g. sub2ind,
# interp1) have been replaced with Julia ones. 

# This software is licensed to C. L. Haley in 2021 under GNU GPL v 2.0

using Interpolations, PCHIPInterpolation

""" 

    sub2ind(A, row, col)

Convert to linear indices 
https://www.mathworks.com/help/matlab/ref/sub2ind.html
"""
function sub2ind(A, row, col) 
    LinearIndices(A)[CartesianIndex(row, col)]
end



"""

    interp1()

1D data interpolation, linear
https://www.mathworks.com/help/matlab/ref/interp1.html?s_tid=doc_ta

# Arguments
- Vector x contains the sample points, and 
- v contains the corresponding values, v(x). 
- Vector xq contains the coordinates of the query points.

# Outputs
- interpolated values of a 1-D function at specific query points using linear interpolation

# Example usage

```
interp1(LinRange(0, 1, 10), log10.(LinRange(0, 1, 10)), [0.5, 0.9], :linear)
```

"""
function interp1(x, v, xq, method = :linear)
    #(method != :linear)||(method != :pchip) && error("Linear Interpolation is the only one implemented here.")
    interpolant = (method == :linear) ? LinearInterpolation(x, v) : PCHIPInterpolation.Interpolator(x, v)
    return interpolant.(xq)
end

"""

    gamini(data, folding)

# Arguments
-`data`: some data vector
-`folding`: The replication factor for every element of the data; if a scalar thsi applies to all of the elements (default 3), if zero or negative, no replication occurs

# Outputs
-`bigger`

# Example usage
```
a, b = degamini(gamini([1, 2, 3, 1, 4, 5], [1, 2, 3, 2, 4, 2]))
```
One gets [1, 2, 2, 3, 3, 3, 1, 1, 4, 4, 4, 4, 5, 5] as the intermediate result.

# See also
@degamini, @gamini2
"""
function gamini(data, folding = 3)
    # Copy a single scalar folding for all elements in the array
    if prod(size(folding)) == 1
        folding = repeat(folding, prod(size(data)), 1)
    end
    # Never had a replication factor of zero before
    # But only if they are of the same length!
    if length(data) == length(folding)
        ind = findall(f->f>0, folding)
        data = data[ind]
        folding = folding[ind]
    end
    # Rearrange
    data = data[:]
    folding = folding[:]
    !(size(data) == size(folding)) && error("Sizes of input and folding must be the same.")
    gelp = zeros(Int64, sum(folding))
    gelp[vcat(1, cumsum(folding[1:end-1], dims = 1) .+ 1)] .= 1
    if !isempty(data)
        bigger  = data[cumsum(gelp, dims = 1)]
    else
        bigger = []
    end
    return bigger
end

""" 

    degamini(v) 

# Arguments
- `v` A vector with repeated entries

# Outputs
- `dv` The same vector with all the repeat sset to 1
- `foldi` A same-dimensional vector with how many repeats there are
- `be` A matrix with begin and end indices into the original vector

# Example
```
dv, foldi, be = degamini([1, 2, 2, 3])
```
answer is dv = [1, 2, 3]; foldi = [1, 2, 1]; be = [1 1; 2 3; 4 4]
"""
function degamini(v)
    v = v[:]
    indi = vcat(1, findall(diff(v) .!= 0) .+ 1) 
    dv = v[indi]
    foldi = diff(vcat(indi, length(v) + 1))
    (length(v) != sum(foldi)) && error("Incorrect")
    be = hcat(indi[:], cumsum(foldi[:]))
    return dv, foldi, be
end

""" 
    
    matranges(ranges)
    
Makes an index vector with monotonically increasing indices between pairs of numbers supplied as input.
From slepian_alpha
    
# Arguments
- `ranges`

    # Outputs
    
# Example
```
matranges([1, 4, 1, 2, -1, 2])
# answer is [1, 2, 3, 4, 1, 2, -1, 0, 1, 2]
```
"""
function matranges(ranges)
    ranges = ranges[:]
    (mod(length(ranges), 2) != 0) && error("Ranges must form pairs")
    lower = ranges[1:2:end]
    upper = ranges[2:2:end]
    hulp1 = ones(sum(upper .- lower .+ 1))
    hulp2 = cumsum(upper[1:end-1] .- lower[1:end-1] .+ 1)
    hulp1[hulp2 .+ 1] = lower[2:end] .- upper[1:(end-1)]   
    hulp1[1] = ranges[1]
    return Int64.(cumsum(hulp1))
end

""" 

    randcirc(xm, ym, r,dr, N)

# Arguments
- xm horizontal positon of the center
- ym vertical position of the center
- `r` radius
- `dr` size of random perturbations around the radius
- `N` number of random spike points

# Outputs
- `x` x-coordinate
- `y` y-coordinate

# Related
@blob

"""
function randcirc(xm = 0.0, ym = 0.0, r = 1.0, dr = 1, N = 10)
    nr = 100
    r = r*ones(N)
    r = r + 2dr*(rand(N) .- 0.5)
    t = LinRange(0, 2*pi*N/(N+1), N)
    
    r = vcat(r,r,r)
    t = vcat(t .- 2 * pi, t, t .+ 2 *pi)
    
    tt = LinRange(0, 2*pi, nr)
    rr = interp1(t, r, collect(tt), :pchip)
    
    x = @. xm + rr * cos(tt)
    y = @. ym + rr * sin(tt)
    
    return x, y
end

""" 
    blob(N, Nj)

Makes moving picture of a random blob by superposition of random circles

# Arguments
- `N` number of loops for, and if movie, default = 100
- `Nj` smoothness, roughly (?)

# Outputs
- `x`
- `y` 
"""
function blob(N = 100, Nj = 10)
    xold, yold = randcirc(0, 0, 1, 1, 10) ### Help: what is randcirc
    r = LinRange(0, 1, Nj+1)
    rm = 1 .- r
    xx = []
    yy = []
    for index = 1:N
        x,y = randcirc(0, 0, 1, 0.2, 10)
        for j = 1:Nj
            push!(xx, x * r[j] + xold * rm[j])
            push!(yy, y * r[j] + yold * rm[j])
        end
        xold = x
        yold = y
    end
    return xx, yy
end

""" 

    phicurve(thph, th)

Adapted from slepian_alpha; finds the longitude crossings and thus integration domains of a closed curve parameterized in colatitude/longitude space at certian query points of colatitude.

# Arguments
- `thph::` Colatitude/longitude of the closed curve (degrees)
- `th::` Colatitude at which crossings are required (degrees)

# Outputs
- `phint` A matrix with crossings/intervals and zeros of dimensions MxN where M = length(th) and N can be anything depending on the oscillations of the curve
- `thp` Colatitude matrix for hatched plotting, if possible
- `php` Longitude matrix for hatched plotting, if possible
- `forreal` Indices of the ones that are real (could be at zero)

Depends on: @sub2ind, @interp1, @degamini, @matranges, and possibly @blob (demo2)
"""
function phicurve(thph, th)
    # Fore every th, find the relevant phint
    xmth = repeat(thph[:,1], 1, length(th)) .- repeat(th', length(thph[:,1]), 1)
    dsx = diff(sign.(xmth), dims = 1)
    # It can be the one before, or after the crossing
    # colf, colj = findall(x -> (x != 0.0), dsx)
    col = findall(x -> (x != 0.0), dsx)
    colf, colj = (map(c -> Tuple(c)[1], col), map(c -> Tuple(c)[2], col))
    # colr = sub2ind(dsx, colf, colj) 
    colr = LinearIndices(dsx)[col]
    # This returns the one on the negative side of hte line
    # add one to colx if the difference is -2; add one to colx2 if the difference is 2
    colx = colf .+ (dsx[colr] .== -2) # DOT was missing here in parentheses
    colx2 = colf .+ (dsx[colr] .== 2)
    L = length(colx)
    (L%2 == 1) && error("Cannot find pairs of crossings.")
    phint = zeros(L)
    # Then one point was exactly hit, this is the thN or the thS case
    if length(colx) == 2 && colx == colx2
        phint = thph[hcat(colx[2], colx2[2]), 2]
        thp = hcat(th, th)
        php = phint
    else
        for ond = 1:L
            # In case you have a node immediately followed by a crossing
            phint[ond] = (colx[ond] == colx2[ond]) ? NaN : interp1(xmth[hcat(colx[ond], colx2[ond]), colj[ond]][:], thph[hcat(colx[ond], colx2[ond]),2][:], 0, :linear) 
        end
    end
    # If the NaNs are not consecutive pairs, get special case
    # Now rearrange back othte number of requested points
    # There could be points with more or less than two crossings
    # Maximum number of times a crossing is repeated
    a, b =  degamini(colj) 
    rowj = copy(colj)
    colj = matranges(Int64.(reshape(hcat(ones(length(b)), b)', length(b)*2, 1)))
    pint = repeat([NaN], length(th), maximum(b))
    subsi = (colj .- 1) .* length(th) .+ rowj # Linear index
    pint[subsi] = phint
    wt, thp = (length(b) == length(th)) ? (0, reshape(gamini(th, b), (2, Int64(length(phint)/2)))) : (1, [])
    # Need to sort since contour may be given in any order
    phint = sort(pint, dims = 2)
    php = (wt == 0) ? reshape(phint[subsi], (2, Int64(length(colj)/2))) : [] # line 95
    # Make them zero so the integral doesn't do anything
    forreal = map(x -> !isnan(x), phint)
    phint[findall(isnan, phint)] .= 0.0
    # note: can use (demo2)
    # x,y=blob(1,1); thph=hcat(y[:],x[:]); Nth = ceil(rand*300); th=LinRange(minimum(thph[:,1]), maximum(thph[:,1]), Nth)
    return phint, thp, php, forreal
end

""" 
    quadpts(qx, Nqx, qy, forreal, xints, wqx, wqy)

Scale the quaduature points and weights 

# Arguments
- `qx` 
- `Nqx`
- `qy`
- `forreal`
- `xints`
- `wqx`
- `wqy`

# Outputs
- `QX` Quadrature points in x
- `QY` Quadrature points in y
- `w` Quadrature weights
- `Nrun` Number of horizontal line segments in the domain

"""
function quadpts(qx, Nqx, qy, forreal, xints, wqx, wqy)
    Nall = Int64(sum(forreal[:])/2)
    # Initialize quadrature points
    QX = repeat([NaN], Nqx, Nall)
    QY = repeat([NaN], Nqx, Nall)
    w = repeat([NaN], Nqx, Nall)
    Nrun = 0
    for yindex = 1:length(qy)
        # How many horizontal intervals at this y?
        Nxint = Int64(sum(forreal[yindex,:])/2)
        abset = reshape(xints[yindex, 1:Nxint*2], 2, Nxint)
        # beginning and endpoints of those intervals
        a = abset[1,:]
        b = abset[2,:]
        # produce the scaled nodes, a column per interval
        exs = repeat(a', length(qx), 1) .+ (qx .+ 1)/2 .* (b .- a)'
        # scatter!(p, exs, -0.5*ones(32), marker = :o)
        # And produce the weights that go with it 
        # and multiply it with the current y node and weight 
        # and produce copies of the y nodes, as many as needed 
        # and put them into two big vectors of weights and points
        QX[:, (Nrun+1):(Nrun+Nxint)]= exs
        QY[:, (Nrun+1):(Nrun+Nxint)] = repeat([qy[yindex]], size(exs)...)
        w[:, (Nrun+1):(Nrun+Nxint)] = wqx.*((b.-a)'/2)*wqy[yindex]
        Nrun = Nrun + Nxint
    end
    return QX, QY, w, Nrun
end

sign(x) = (x >= 0.0) ? 1 : -1

function get_quadrature_nodes_2D(x, y, Nqx = 32, Nqy = 32)
    #println("Gauss-Legendre method with ($Nqx, $Nqy) integration nodes.")

    # Calculate blanks, xaxis GL points regardless of the data range
    qx, wqx = FastGaussQuadrature.gausslegendre(Nqx)
    qy, wqy = FastGaussQuadrature.gausslegendre(Nqy)
    
    mn, mx = (minimum(y), maximum(y))
    # Scale the y-quadrature points to the min-max interval
    th = qy*(abs(mx - mn)/2) .+ (mx + mn)/2

    # For each of these points in y find the appropriate ranges of x's
    xints,yp,xp,forreal = phicurve(hcat(y, x), th) 
    QX, QY, w, Nrun = quadpts(qx, Nqx, th, forreal, xints, wqx, wqy)

    return QX, QY, w, Nrun
end

#  GNU GENERAL PUBLIC LICENSE
                       # Version 2, June 1991

 # Copyright (C) 1989, 1991 Free Software Foundation, Inc.,
 # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 # Everyone is permitted to copy and distribute verbatim copies
 # of this license document, but changing it is not allowed.

                            # Preamble

  # The licenses for most software are designed to take away your
# freedom to share and change it.  By contrast, the GNU General Public
# License is intended to guarantee your freedom to share and change free
# software--to make sure the software is free for all its users.  This
# General Public License applies to most of the Free Software
# Foundation's software and to any other program whose authors commit to
# using it.  (Some other Free Software Foundation software is covered by
# the GNU Lesser General Public License instead.)  You can apply it to
# your programs, too.

  # When we speak of free software, we are referring to freedom, not
# price.  Our General Public Licenses are designed to make sure that you
# have the freedom to distribute copies of free software (and charge for
# this service if you wish), that you receive source code or can get it
# if you want it, that you can change the software or use pieces of it
# in new free programs; and that you know you can do these things.

  # To protect your rights, we need to make restrictions that forbid
# anyone to deny you these rights or to ask you to surrender the rights.
# These restrictions translate to certain responsibilities for you if you
# distribute copies of the software, or if you modify it.

#   For example, if you distribute copies of such a program, whether
# gratis or for a fee, you must give the recipients all the rights that
# you have.  You must make sure that they, too, receive or can get the
# source code.  And you must show them these terms so they know their
# rights.

  # We protect your rights with two steps: (1) copyright the software, and
# (2) offer you this license which gives you legal permission to copy,
# distribute and/or modify the software.

  # Also, for each author's protection and ours, we want to make certain
# that everyone understands that there is no warranty for this free
# software.  If the software is modified by someone else and passed on, we
# want its recipients to know that what they have is not the original, so
# that any problems introduced by others will not reflect on the original
# authors' reputations.

  # Finally, any free program is threatened constantly by software
# patents.  We wish to avoid the danger that redistributors of a free
# program will individually obtain patent licenses, in effect making the
# program proprietary.  To prevent this, we have made it clear that any
# patent must be licensed for everyone's free use or not licensed at all.

  # The precise terms and conditions for copying, distribution and
# modification follow.

                    # GNU GENERAL PUBLIC LICENSE
   # TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

  # 0. This License applies to any program or other work which contains
# a notice placed by the copyright holder saying it may be distributed
# under the terms of this General Public License.  The "Program", below,
# refers to any such program or work, and a "work based on the Program"
# means either the Program or any derivative work under copyright law:
# that is to say, a work containing the Program or a portion of it,
# either verbatim or with modifications and/or translated into another
# language.  (Hereinafter, translation is included without limitation in
# the term "modification".)  Each licensee is addressed as "you".

# Activities other than copying, distribution and modification are not
# covered by this License; they are outside its scope.  The act of
# running the Program is not restricted, and the output from the Program
# is covered only if its contents constitute a work based on the
# Program (independent of having been made by running the Program).
# Whether that is true depends on what the Program does.

  # 1. You may copy and distribute verbatim copies of the Program's
# source code as you receive it, in any medium, provided that you
# conspicuously and appropriately publish on each copy an appropriate
# copyright notice and disclaimer of warranty; keep intact all the
# notices that refer to this License and to the absence of any warranty;
# and give any other recipients of the Program a copy of this License
# along with the Program.

# You may charge a fee for the physical act of transferring a copy, and
# you may at your option offer warranty protection in exchange for a fee.

  # 2. You may modify your copy or copies of the Program or any portion
# of it, thus forming a work based on the Program, and copy and
# distribute such modifications or work under the terms of Section 1
# above, provided that you also meet all of these conditions:

    # a) You must cause the modified files to carry prominent notices
    # stating that you changed the files and the date of any change.

    # b) You must cause any work that you distribute or publish, that in
    # whole or in part contains or is derived from the Program or any
    # part thereof, to be licensed as a whole at no charge to all third
    # parties under the terms of this License.

    # c) If the modified program normally reads commands interactively
    # when run, you must cause it, when started running for such
    # interactive use in the most ordinary way, to print or display an
    # announcement including an appropriate copyright notice and a
    # notice that there is no warranty (or else, saying that you provide
    # a warranty) and that users may redistribute the program under
    # these conditions, and telling the user how to view a copy of this
    # License.  (Exception: if the Program itself is interactive but
    # does not normally print such an announcement, your work based on
    # the Program is not required to print an announcement.)

# These requirements apply to the modified work as a whole.  If
# identifiable sections of that work are not derived from the Program,
# and can be reasonably considered independent and separate works in
# themselves, then this License, and its terms, do not apply to those
# sections when you distribute them as separate works.  But when you
# distribute the same sections as part of a whole which is a work based
# on the Program, the distribution of the whole must be on the terms of
# this License, whose permissions for other licensees extend to the
# entire whole, and thus to each and every part regardless of who wrote it.

# Thus, it is not the intent of this section to claim rights or contest
# your rights to work written entirely by you; rather, the intent is to
# exercise the right to control the distribution of derivative or
# collective works based on the Program.

# In addition, mere aggregation of another work not based on the Program
# with the Program (or with a work based on the Program) on a volume of
# a storage or distribution medium does not bring the other work under
# the scope of this License.

  # 3. You may copy and distribute the Program (or a work based on it,
# under Section 2) in object code or executable form under the terms of
# Sections 1 and 2 above provided that you also do one of the following:

    # a) Accompany it with the complete corresponding machine-readable
    # source code, which must be distributed under the terms of Sections
    # 1 and 2 above on a medium customarily used for software interchange; or,

    # b) Accompany it with a written offer, valid for at least three
    # years, to give any third party, for a charge no more than your
    # cost of physically performing source distribution, a complete
    # machine-readable copy of the corresponding source code, to be
    # distributed under the terms of Sections 1 and 2 above on a medium
    # customarily used for software interchange; or,

    # c) Accompany it with the information you received as to the offer
    # to distribute corresponding source code.  (This alternative is
    # allowed only for noncommercial distribution and only if you
    # received the program in object code or executable form with such
    # an offer, in accord with Subsection b above.)

# The source code for a work means the preferred form of the work for
# making modifications to it.  For an executable work, complete source
# code means all the source code for all modules it contains, plus any
# associated interface definition files, plus the scripts used to
# control compilation and installation of the executable.  However, as a
# special exception, the source code distributed need not include
# anything that is normally distributed (in either source or binary
# form) with the major components (compiler, kernel, and so on) of the
# operating system on which the executable runs, unless that component
# itself accompanies the executable.

# If distribution of executable or object code is made by offering
# access to copy from a designated place, then offering equivalent
# access to copy the source code from the same place counts as
# distribution of the source code, even though third parties are not
# compelled to copy the source along with the object code.

  # 4. You may not copy, modify, sublicense, or distribute the Program
# except as expressly provided under this License.  Any attempt
# otherwise to copy, modify, sublicense or distribute the Program is
# void, and will automatically terminate your rights under this License.
# However, parties who have received copies, or rights, from you under
# this License will not have their licenses terminated so long as such
# parties remain in full compliance.

#   5. You are not required to accept this License, since you have not
# signed it.  However, nothing else grants you permission to modify or
# distribute the Program or its derivative works.  These actions are
# prohibited by law if you do not accept this License.  Therefore, by
# modifying or distributing the Program (or any work based on the
# Program), you indicate your acceptance of this License to do so, and
# all its terms and conditions for copying, distributing or modifying
# the Program or works based on it.

  # 6. Each time you redistribute the Program (or any work based on the
# Program), the recipient automatically receives a license from the
# original licensor to copy, distribute or modify the Program subject to
# these terms and conditions.  You may not impose any further
# restrictions on the recipients' exercise of the rights granted herein.
# You are not responsible for enforcing compliance by third parties to
# this License.

  # 7. If, as a consequence of a court judgment or allegation of patent
# infringement or for any other reason (not limited to patent issues),
# conditions are imposed on you (whether by court order, agreement or
# otherwise) that contradict the conditions of this License, they do not
# excuse you from the conditions of this License.  If you cannot
# distribute so as to satisfy simultaneously your obligations under this
# License and any other pertinent obligations, then as a consequence you
# may not distribute the Program at all.  For example, if a patent
# license would not permit royalty-free redistribution of the Program by
# all those who receive copies directly or indirectly through you, then
# the only way you could satisfy both it and this License would be to
# refrain entirely from distribution of the Program.

# If any portion of this section is held invalid or unenforceable under
# any particular circumstance, the balance of the section is intended to
# apply and the section as a whole is intended to apply in other
# circumstances.

# It is not the purpose of this section to induce you to infringe any
# patents or other property right claims or to contest validity of any
# such claims; this section has the sole purpose of protecting the
# integrity of the free software distribution system, which is
# implemented by public license practices.  Many people have made
# generous contributions to the wide range of software distributed
# through that system in reliance on consistent application of that
# system; it is up to the author/donor to decide if he or she is willing
# to distribute software through any other system and a licensee cannot
# impose that choice.

# This section is intended to make thoroughly clear what is believed to
# be a consequence of the rest of this License.

#   8. If the distribution and/or use of the Program is restricted in
# certain countries either by patents or by copyrighted interfaces, the
# original copyright holder who places the Program under this License
# may add an explicit geographical distribution limitation excluding
# those countries, so that distribution is permitted only in or among
# countries not thus excluded.  In such case, this License incorporates
# the limitation as if written in the body of this License.

  # 9. The Free Software Foundation may publish revised and/or new versions
# of the General Public License from time to time.  Such new versions will
# be similar in spirit to the present version, but may differ in detail to
# address new problems or concerns.

# Each version is given a distinguishing version number.  If the Program
# specifies a version number of this License which applies to it and "any
# later version", you have the option of following the terms and conditions
# either of that version or of any later version published by the Free
# Software Foundation.  If the Program does not specify a version number of
# this License, you may choose any version ever published by the Free Software
# Foundation.

  # 10. If you wish to incorporate parts of the Program into other free
# programs whose distribution conditions are different, write to the author
# to ask for permission.  For software which is copyrighted by the Free
# Software Foundation, write to the Free Software Foundation; we sometimes
# make exceptions for this.  Our decision will be guided by the two goals
# of preserving the free status of all derivatives of our free software and
# of promoting the sharing and reuse of software generally.

                            # NO WARRANTY

  # 11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
# FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
# OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
# PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
# OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
# TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
# PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
# REPAIR OR CORRECTION.

  # 12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
# WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
# REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
# INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
# OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
# TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
# YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
# PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGES.

                     # END OF TERMS AND CONDITIONS

            # How to Apply These Terms to Your New Programs

#   If you develop a new program, and you want it to be of the greatest
# possible use to the public, the best way to achieve this is to make it
# free software which everyone can redistribute and change under these terms.

  # To do so, attach the following notices to the program.  It is safest
# to attach them to the start of each source file to most effectively
# convey the exclusion of warranty; and each file should have at least
# the "copyright" line and a pointer to where the full notice is found.

    # <one line to give the program's name and a brief idea of what it does.>
    # Copyright (C) <year>  <name of author>

    # This program is free software; you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation; either version 2 of the License, or
    # (at your option) any later version.

    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.

    # You should have received a copy of the GNU General Public License along
    # with this program; if not, write to the Free Software Foundation, Inc.,
    # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# Also add information on how to contact you by electronic and paper mail.

# If the program is interactive, make it output a short notice like this
# when it starts in an interactive mode:

    # Gnomovision version 69, Copyright (C) year name of author
    # Gnomovision comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    # This is free software, and you are welcome to redistribute it
    # under certain conditions; type `show c' for details.

# The hypothetical commands `show w' and `show c' should show the appropriate
# parts of the General Public License.  Of course, the commands you use may
# be called something other than `show w' and `show c'; they could even be
# mouse-clicks or menu items--whatever suits your program.

# You should also get your employer (if you work as a programmer) or your
# school, if any, to sign a "copyright disclaimer" for the program, if
# necessary.  Here is a sample; alter the names:

  # Yoyodyne, Inc., hereby disclaims all copyright interest in the program
  # `Gnomovision' (which makes passes at compilers) written by James Hacker.

  # <signature of Ty Coon>, 1 April 1989
  # Ty Coon, President of Vice

# This General Public License does not permit incorporating your program into
# proprietary programs.  If your program is a subroutine library, you may
# consider it more useful to permit linking proprietary applications with the
# library.  If this is what you want to do, use the GNU Lesser General
# Public License instead of this License.
